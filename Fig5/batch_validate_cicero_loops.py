#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
batch validate Cicero Co-accessibility Pairs and Hi-C Loops
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
from datetime import datetime

SAMPLES = [ "Endothelial cell_KO",'Atrial cardiomyocyte_KO','Immune cell_KO','Macrophage_KO']
CHROMOSOMES = [f"chr{i}" for i in range(1, 20) if i != 4]  
SELECTED_SLOP = 20000  # 20kb

BASE_DIR = Path("/media/ggj/ggj/CJY/nature_WXY/Hi-c/analysis/dots")
# CICERO_BASE = Path("/media/ggj/ggj/CJY/cuttag/scATAC/tissue/Heart/cicero_downsample/result/process")
CICERO_BASE = Path("/media/ggj/ggj/CJY/cuttag/scATAC/tissue/Heart/cicero_celltype/result/process")
LOG_FILE = BASE_DIR / "result" / "batch_validation_log.txt"

def log_message(message, log_file=LOG_FILE):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_msg = f"[{timestamp}] {message}"
    print(log_msg)
    with open(log_file, "a") as f:
        f.write(log_msg + "\n")

def is_valid_in_interaction_region(c1, c2, interaction_regions_df, slop=0):
    """check if cicero pair's two centers are in loops"""
    if c1 > c2:
        c1, c2 = c2, c1
    
    for idx, region in interaction_regions_df.iterrows():
        region_start = region['start'] - slop
        region_end = region['end'] + slop
        
        if (region_start <= c1 <= region_end) and (region_start <= c2 <= region_end):
            return True, region['type'], idx
    
    return False, None, -1

def is_valid_pair_with_slop(c1, c2, loops_df, slop=0):
    """check if one center is in loops"""
    if c1 > c2:
        c1, c2 = c2, c1
    
    loops = loops_df[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']].copy()
    loops.columns = ['chr1', 's1', 'e1', 'chr2', 's2', 'e2']
    
    for idx, loop in loops.iterrows():
        s1_ext = loop['s1'] - slop
        e1_ext = loop['e1'] + slop
        s2_ext = loop['s2'] - slop
        e2_ext = loop['e2'] + slop
        
        if (s1_ext <= c1 <= e1_ext) and (s2_ext <= c2 <= e2_ext):
            return True, 'single_loop', -1
        if (s2_ext <= c1 <= e2_ext) and (s1_ext <= c2 <= e1_ext):
            return True, 'single_loop', -1
    return False, None, -1

def process_chromosome(sample, chrom):
    """processing each chr"""

    sample_suffix = sample.split("_", 1)[1] if "_" in sample else sample
    cicero_file = CICERO_BASE / f"{sample}.txt"
    dots_file = BASE_DIR / f"{sample_suffix}_10bin" / f"{sample_suffix}_noXY_10kb_{chrom}_dots.bedpe"
    output_dir = BASE_DIR / "result" / sample
    output_file = output_dir / f"cicero_hic_validation_{sample}_{chrom}_interaction_regions.tsv"
    valid_output = output_dir / f"cicero_hic_validation_{sample}_{chrom}_interaction_regions_valid_only.tsv"
    regions_output = output_dir / f"interaction_regions_{sample}_{chrom}.tsv"

    if not cicero_file.exists():
        log_message(f"  error: Cicero 文件不存在: {cicero_file}")
        return False
    
    if not dots_file.exists():
        log_message(f"  error: Dots 文件不存在: {dots_file}")
        return False

    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        log_message(f"  load Hi-C dots: {dots_file}")
        dots = pd.read_csv(dots_file, sep="\t", header=0)
        dots.columns = dots.columns.str.strip()
        
        log_message(f"  load Cicero pairs: {cicero_file}")
        cicero_cols = [
            'id', 'peak1', 'peak2', 'coaccess', 'is_promoter1', 'is_promoter2', 
            'interaction_count', 'gene1', 'gene2', 'unique_id', 'center1', 'center2',
            'distance', 'dist_class', 'annot1_detail', 'annot1_type', 'annot1_gene',
            'annot2_detail', 'annot2_type', 'annot2_gene', 'pair_type'
        ]
        cicero = pd.read_csv(cicero_file, sep="\t", header=None, names=cicero_cols,index_col=False)
        cicero['chr'] = cicero['peak1'].str.split('_').str[0]
        cicero_chr = cicero[cicero['chr'] == chrom].copy()
        cicero_chr['center1'] = pd.to_numeric(cicero_chr['center1'], errors='coerce')
        cicero_chr['center2'] = pd.to_numeric(cicero_chr['center2'], errors='coerce')
        cicero_chr = cicero_chr.dropna(subset=['center1', 'center2'])
        
        log_message(f"  Hi-C loops: {len(dots):,}, Cicero pairs ({chrom}): {len(cicero_chr):,}")
        
        if len(dots) == 0:
            log_message(f"  warning: {chrom} no Hi-C loops")
            return False
        
        if len(cicero_chr) == 0:
            log_message(f"  warning: {chrom} no Cicero pairs")
            return False

        valid_loops = dots[
            (dots['la_exp.donut.qval'] < 0.1) &
            (dots['la_exp.donut.value'] > 5) &
            (dots['count'] > 10)
        ].copy()
        
        if len(valid_loops) == 0:
            log_message(f"  warnig: {chrom} no high quality loops")
            return False
        
        valid_loops['anchor1_center'] = (valid_loops['start1'] + valid_loops['end1']) / 2
        valid_loops['anchor2_center'] = (valid_loops['start2'] + valid_loops['end2']) / 2
        
        log_message(f"  high quality loops: {len(valid_loops):,}/{len(dots):,}")

        proximity_threshold = 50000  # 50kb
        neighboring_pairs = []
        
        for i, loop1 in valid_loops.iterrows():
            for j, loop2 in valid_loops.iterrows():
                if i >= j:
                    continue
                
                dists = [
                    abs(loop1['anchor1_center'] - loop2['anchor1_center']),
                    abs(loop1['anchor1_center'] - loop2['anchor2_center']),
                    abs(loop1['anchor2_center'] - loop2['anchor1_center']),
                    abs(loop1['anchor2_center'] - loop2['anchor2_center'])
                ]
                min_dist = min(dists)
                
                if min_dist <= proximity_threshold:
                    min_idx = dists.index(min_dist)
                    if min_idx == 0:
                        a1_start, a1_end = loop1['start1'], loop1['end1']
                        a2_start, a2_end = loop2['start1'], loop2['end1']
                    elif min_idx == 1:
                        a1_start, a1_end = loop1['start1'], loop1['end1']
                        a2_start, a2_end = loop2['start2'], loop2['end2']
                    elif min_idx == 2:
                        a1_start, a1_end = loop1['start2'], loop1['end2']
                        a2_start, a2_end = loop2['start1'], loop2['end1']
                    else:
                        a1_start, a1_end = loop1['start2'], loop1['end2']
                        a2_start, a2_end = loop2['start2'], loop2['end2']
                    
                    neighboring_pairs.append({
                        'chr': loop1['chrom1'],
                        'start': min(a1_start, a2_start) - 20000,
                        'end': max(a1_end, a2_end) + 20000,
                        'type': 'neighboring_loops',
                        'num_loops': 2
                    })
        
        all_anchors = []
        for idx, row in valid_loops.iterrows():
            all_anchors.append({
                'chr': row['chrom1'],
                'start': row['start1'],
                'end': row['end1'],
                'center': row['anchor1_center'],
                'loop_idx': idx
            })
            all_anchors.append({
                'chr': row['chrom2'],
                'start': row['start2'],
                'end': row['end2'],
                'center': row['anchor2_center'],
                'loop_idx': idx
            })
        
        anchors_df = pd.DataFrame(all_anchors)
        anchor_overlap_threshold = 10000  # 10kb
        hub_regions = []
        processed_anchors = set()
        
        for anchor in anchors_df.itertuples():
            anchor_key = (anchor.chr, anchor.center)
            if anchor_key in processed_anchors:
                continue
            
            overlapping = anchors_df[
                (anchors_df['chr'] == anchor.chr) &
                (abs(anchors_df['center'] - anchor.center) <= anchor_overlap_threshold) &
                (anchors_df['loop_idx'] != anchor.loop_idx)
            ]
            
            if len(overlapping) > 0:
                all_starts = [anchor.start] + overlapping['start'].tolist()
                all_ends = [anchor.end] + overlapping['end'].tolist()
                hub_regions.append({
                    'chr': anchor.chr,
                    'start': min(all_starts),
                    'end': max(all_ends),
                    'type': 'hub',
                    'num_loops': len(overlapping) + 1
                })
                processed_anchors.add(anchor_key)

        cluster_regions = []
        if 'c_label' in valid_loops.columns:
            for cluster_label in valid_loops['c_label'].unique():
                cluster_loops = valid_loops[valid_loops['c_label'] == cluster_label]
                
                if len(cluster_loops) > 1:  
                    all_starts = pd.concat([cluster_loops['start1'], cluster_loops['start2']])
                    all_ends = pd.concat([cluster_loops['end1'], cluster_loops['end2']])
                    
                    cluster_regions.append({
                        'chr': cluster_loops.iloc[0]['chrom1'],
                        'start': all_starts.min(),
                        'end': all_ends.max(),
                        'type': 'cluster',
                        'num_loops': len(cluster_loops)
                    })

        interaction_regions = neighboring_pairs + hub_regions + cluster_regions
        if len(interaction_regions) > 0:
            interaction_regions_df = pd.DataFrame(interaction_regions)
            interaction_regions_df = interaction_regions_df.sort_values(['chr', 'start'])
            interaction_regions_chr = interaction_regions_df[
                interaction_regions_df['chr'] == chrom
            ].copy()
        else:
            interaction_regions_chr = pd.DataFrame(columns=['chr', 'start', 'end', 'type', 'num_loops'])
        
        log_message(f"  interaction region: {len(interaction_regions_chr):,} (neighboring: {len(neighboring_pairs)}, hub: {len(hub_regions)}, cluster: {len(cluster_regions)})")
        
        results = []
        valid_count = 0
        total = len(cicero_chr)
        
        for i, (idx, row) in enumerate(cicero_chr.iterrows()):
            c1, c2 = row['center1'], row['center2']
            
            if c1 > c2:
                c1, c2 = c2, c1
            
            if len(interaction_regions_chr) > 0:
                is_valid, region_type, matched_idx = is_valid_in_interaction_region(
                    c1, c2, interaction_regions_chr, slop=SELECTED_SLOP
                )
            else:
                is_valid = False
                region_type = None
                matched_idx = -1

            if not is_valid:
                is_valid, region_type, matched_idx = is_valid_pair_with_slop(c1, c2, dots, slop=SELECTED_SLOP)

            if not is_valid:
                region_type = 'none'
            
            results.append({
                'cicero_id': row.get('id', idx),
                'peak1': row['peak1'],
                'peak2': row['peak2'],
                'center1': row['center1'],
                'center2': row['center2'],
                'coaccess': row.get('coaccess', np.nan),
                'distance': row.get('distance', np.nan),
                'pair_type': row.get('pair_type', ''),
                'gene1': row.get('gene1', ''),
                'gene2': row.get('gene2', ''),
                'is_valid': is_valid,
                'matched_region_type': region_type if region_type else 'none',
                'matched_region_idx': matched_idx
            })
            
            if is_valid:
                valid_count += 1
            
            if (i + 1) % 2000 == 0:
                log_message(f"    processing: {i+1}/{total} ({(i+1)/total*100:.1f}%), valid: {valid_count}")
        
        results_df = pd.DataFrame(results)
        valid_rate = valid_count / total * 100 if total > 0 else 0

        results_df.to_csv(output_file, sep="\t", index=False)

        valid_df = results_df[results_df['is_valid']].copy()
        if len(valid_df) > 0:
            valid_df.to_csv(valid_output, sep="\t", index=False)
            log_message(f"  Valid pairs saved: {valid_output} ({len(valid_df)} pairs)")
        
        if len(interaction_regions_chr) > 0:
            interaction_regions_chr.to_csv(regions_output, sep="\t", index=False)
            log_message(f"  saved: {regions_output}")
        
        return True
        
    except Exception as e:
        log_message(f"  error: processing {sample}, {chrom} error: {str(e)}")
        import traceback
        log_message(f"  error: {traceback.format_exc()}")
        return False

def main():
    """main"""
    log_message("="*80)
    log_message("validate Cicero Co-accessibility Pairs and Hi-C Loops")
    log_message(f"samples: {SAMPLES}")
    log_message(f"chr: {CHROMOSOMES}")
    log_message(f"Slop : {SELECTED_SLOP/1000}kb")
    log_message("="*80)
    
    total_tasks = len(SAMPLES) * len(CHROMOSOMES)
    completed = 0
    failed = 0
    
    for sample in SAMPLES:
        log_message(f"\nprocessing: {sample}")
        for chrom in CHROMOSOMES:
            success = process_chromosome(sample, chrom)
            completed += 1
            if success:
                log_message(f"  done: {sample}, {chrom} ({completed}/{total_tasks})")
            else:
                failed += 1
                log_message(f"  failed: {sample}, {chrom} ({completed}/{total_tasks})")
    
    log_message("\n" + "="*80)
    log_message(f"all: {total_tasks}")
    log_message(f"done: {total_tasks - failed}")
    log_message(f"failed: {failed}")
    log_message("="*80)

if __name__ == "__main__":
    main()
