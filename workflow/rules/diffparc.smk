#added some hints for eventual bids-derivatives naming (e.g. space, label, type(dseg, mask, probseg)..)

#space-T1w (native), dseg
rule combine_lr_hcp:
    input:
        lh = bids(root='results/hcp_mmp',subject='{subject}',hemi='L',label='hcpmmp',space='native',suffix='dseg.nii.gz'),
        rh = bids(root='results/hcp_mmp',subject='{subject}',hemi='R',label='hcpmmp',space='native',suffix='dseg.nii.gz'),
    output:
        lh_rh = bids(root='results/diffparc',subject='{subject}',space='individual',label=config['targets_atlas_name'],suffix='dseg.nii.gz')
    container: config['singularity']['neuroglia']
    log: 'logs/combine_lr_hcp/{subject}.log'
    group: 'participant1'
    shell: 'fslmaths {input.lh} -max {input.rh} {output.lh_rh} &> {log}'


#space-{template}, probseg
rule get_template_seed:
    input: 
        seed = lambda wildcards: config['template_prob_seg'][wildcards.seed],
    output: bids(root='results/diffparc',template='{template}',label='{seed}',suffix='probseg.nii.gz')
    container: config['singularity']['neuroglia']
    log: 'logs/get_template_seed/{template}_{seed}.log'
    group: 'group0'
    shell:
        'cp -v {input} {output} &> {log}'
 

# mask
rule binarize_template_seed:
    input: 
        seed = bids(root='results/diffparc',template='{template}',label='{seed}',suffix='probseg.nii.gz')
    params:
        threshold = config['prob_seg_threshold']
    output: 
        mask = bids(root='results/diffparc',template='{template}',label='{seed}',suffix='mask.nii.gz')
    container: config['singularity']['neuroglia']
    log: 'logs/binarize_template_seed/{template}_{seed}.log'
    group: 'group0'
    shell:
        'fslmaths {input} -thr {params.threshold} {output} &> {log}'
        




#transform probabilistic seed to subject
#space-T1w,  probseg
rule transform_to_subject:
    input: 
        seed = bids(root='results/diffparc',template='{template}',label='{seed}',suffix='probseg.nii.gz'),
        affine =  config['ants_affine_mat'],
        invwarp =  config['ants_invwarp_nii'],
        ref = bids(root='results/diffparc',subject='{subject}',space='individual',label=config['targets_atlas_name'],suffix='dseg.nii.gz')
    output: 
        seed = bids(root='results/diffparc',subject='{subject}',space='individual',label='{seed}',from_='{template}',suffix='probseg.nii.gz'),
    envmodules: 'ants'
    container: config['singularity']['neuroglia']
    log: 'logs/transform_to_subject/{template}_sub-{subject}_{seed}.log'
    group: 'participant1'
    threads: 8
    shell:
        'antsApplyTransforms -d 3 --interpolation Linear -i {input.seed} -o {output} -r {input.ref} -t [{input.affine},1] -t {input.invwarp} &> {log}'

#create brainmask from bedpost data, and resample to chosen resolution
#space-T1w res-? mask
rule resample_brainmask:
    input: 
        dwi = join(config['fsl_bedpost_dir'],config['bedpost_mean_s0_samples']),
    params:
        seed_resolution = config['probtrack']['seed_resolution']
    output:
        mask = bids(root='results/diffparc',subject='{subject}',label='brain',suffix='mask.nii.gz'),
        mask_res = bids(root='results/diffparc',subject='{subject}',label='brain',res='dwi',suffix='mask.nii.gz'),
    container: config['singularity']['neuroglia']
    log: 'logs/resample_brainmask/sub-{subject}.log'
    group: 'participant1'
    shell:
        'fslmaths {input.dwi} -bin {output.mask} &&'
        'mri_convert {output.mask} -vs {params.seed_resolution} {params.seed_resolution} {params.seed_resolution} {output.mask_res} -rt nearest &> {log}'



#resample target segs to the match resampled dwi mask
#space-T1w res-? dseg
rule resample_targets:
    input: 
        mask_res = bids(root='results/diffparc',subject='{subject}',label='brain',res='dwi',suffix='mask.nii.gz'),
        targets = bids(root='results/diffparc',subject='{subject}',space='individual',label=config['targets_atlas_name'],suffix='dseg.nii.gz')
    output:
        targets_res = bids(root='results/diffparc',subject='{subject}',space='individual',label=config['targets_atlas_name'],res='dwi',suffix='dseg.nii.gz'),
    container: config['singularity']['neuroglia']
    log: 'logs/resample_targets/sub-{subject}.log'
    group: 'participant1'
    shell:
        'reg_resample -flo {input.targets} -res {output.targets_res} -ref {input.mask_res} -NN 0  &> {log}'



#resamples seed seg to match resampled dwi mask
#space-T1w res=? probseg
rule resample_seed:
    input: 
        mask_res = bids(root='results/diffparc',subject='{subject}',label='brain',res='dwi',suffix='mask.nii.gz'),
        seed = bids(root='results/diffparc',subject='{subject}',space='individual',label='{seed}',from_='{template}',suffix='probseg.nii.gz')
    output:
        seed_res = bids(root='results/diffparc',subject='{subject}',space='individual',label='{seed}',from_='{template}',res='dwi',suffix='probseg.nii.gz')
    container: config['singularity']['neuroglia']
    log: 'logs/resample_seed/{template}_sub-{subject}_{seed}.log'
    group: 'participant1'
    shell:
        #linear interp here now, since probabilistic seg
        'reg_resample -flo {input.seed} -res {output.seed_res} -ref {input.mask_res}  &> {log}'

# space-T1w mask
rule binarize_subject_seed:
    input: 
        seed_res = bids(root='results/diffparc',subject='{subject}',space='individual',label='{seed}',from_='{template}',res='dwi',suffix='probseg.nii.gz')
    params:
        threshold = config['prob_seg_threshold']
    output: 
        seed_thr = bids(root='results/diffparc',subject='{subject}',space='individual',label='{seed}',from_='{template}',res='dwi',suffix='mask.nii.gz')
    container: config['singularity']['neuroglia']
    log: 'logs/binarize_subject_seed/{template}_sub-{subject}_{seed}.log'
    container: config['singularity']['neuroglia']
    group: 'participant1'
    shell:
        'fslmaths {input} -thr {params.threshold} {output} &> {log}'
        


#space-T1w, mask 
rule split_targets:
    input: 
        targets = bids(root='results/diffparc',subject='{subject}',space='individual',label=config['targets_atlas_name'],res='dwi',suffix='dseg.nii.gz')
    params:
        target_nums = lambda wildcards: [str(i) for i in range(len(targets))],
        target_seg = lambda wildcards, output: expand('{target_seg_dir}/sub-{subject}_label-{target}_mask.nii.gz',target_seg_dir=output.target_seg_dir,subject=wildcards.subject,target=targets)
    output:
        target_seg_dir = directory(bids(root='results/diffparc',subject='{subject}',suffix='targets'))
    container: config['singularity']['neuroglia']
    log: 'logs/split_targets/sub-{subject}.log'
    threads: 32 
    group: 'participant1'
    shell:  #TODO: could do this in c3d with less effort.. 
        'mkdir -p {output} && parallel  --jobs {threads} fslmaths {input.targets} -thr {{1}} -uthr {{1}} -bin {{2}} &> {log} ::: {params.target_nums} :::+ {params.target_seg}'

#txt
rule gen_targets_txt:
    input:
        target_seg_dir = bids(root='results/diffparc',subject='{subject}',suffix='targets')
    params:
        target_seg = lambda wildcards, input: expand('{target_seg_dir}/sub-{subject}_label-{target}_mask.nii.gz',target_seg_dir=input.target_seg_dir,subject=wildcards.subject,target=targets)
    output:
        target_txt = 'results/diffparc/sub-{subject}/targets.txt'
    log: 'logs/get_targets_txt/sub-{subject}.log'
    group: 'participant1'
    run:
        f = open(output.target_txt,'w')
        for s in params.target_seg:
            f.write(f'{s}\n')
        f.close()

#probtrack dir out
rule run_probtrack:
    input:
        seed_res = bids(root='results/diffparc',subject='{subject}',space='individual',label='{seed}',from_='{template}',res='dwi',suffix='mask.nii.gz'),
        target_txt = rules.gen_targets_txt.output,
        mask = bids(root='results/diffparc',subject='{subject}',label='brain',suffix='mask.nii.gz'),
        target_seg_dir = bids(root='results/diffparc',subject='{subject}',suffix='targets')
    params:
        bedpost_merged = join(config['fsl_bedpost_dir'],config['bedpost_merged_prefix']),
        probtrack_opts = config['probtrack']['opts'],
        out_target_seg = lambda wildcards, output: expand(bids(root=output.probtrack_dir,include_subject_dir=False,prefix='seeds_to',label='{target}',suffix='mask.nii.gz'), target=targets),
        nsamples = config['probtrack']['nsamples'],
        container = '/project/6050199/akhanf/singularity/bids-apps/fsl_6.0.3_cuda9.1.sif'
    output:
        probtrack_dir = directory(bids(root='results/diffparc',subject='{subject}',label='{seed}',from_='{template}',suffix='probtrack'))
    threads: 8
    resources:
        mem_mb = 32000,
        time = 30, #30 mins
        gpus = 1 #1 gpu
    log: 'logs/run_probtrack/sub-{subject}_{seed}_{template}.log'
    group: 'participant1'
    shell:
        'mkdir -p {output.probtrack_dir} && singularity exec -e --nv {params.container} '
        'probtrackx2_gpu --samples={params.bedpost_merged}  --mask={input.mask} --seed={input.seed_res} '
        '--targetmasks={input.target_txt} --seedref={input.seed_res} --nsamples={params.nsamples} '
        '--dir={output.probtrack_dir} {params.probtrack_opts} -V 2  &> {log}'

#check bids-deriv dwi draft (not in main bids yet)
#space-{template}
rule transform_conn_to_template:
    input:
        probtrack_dir = bids(root='results/diffparc',subject='{subject}',label='{seed}',from_='{template}',suffix='probtrack'),
        affine =  config['ants_affine_mat'],
        warp =  config['ants_warp_nii'],
        ref = config['ants_ref_nii']
    params:
        in_connmap_3d = lambda wildcards, input: expand(bids(root=input.probtrack_dir,include_subject_dir=False,subject=wildcards.subject,prefix='seeds_to',label='{target}',suffix='mask.nii.gz'), target=targets),
        out_connmap_3d = lambda wildcards, output: expand(bids(root=output.probtrack_dir,include_subject_dir=False,subject=wildcards.subject,prefix='seeds_to',label='{target}',suffix='mask.nii.gz'), target=targets),
    output:
        probtrack_dir = directory(bids(root='results/diffparc',subject='{subject}',label='{seed}',space='{template}',suffix='probtrack'))
    envmodules: 'ants'
    container: config['singularity']['neuroglia']
    threads: 32
    resources:
        mem_mb = 128000
    log: 'logs/transform_conn_to_template/sub-{subject}_{seed}_{template}.log'
    group: 'participant1'
    shell:
        'mkdir -p {output} && ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1 parallel  --jobs {threads} antsApplyTransforms -d 3 --interpolation Linear -i {{1}} -o {{2}}  -r {input.ref} -t {input.warp} -t {input.affine} &> {log} :::  {params.in_connmap_3d} :::+ {params.out_connmap_3d}' 


#check bids-deriv -- connectivity?
#space-{template}
rule save_connmap_template_npz:
    input:
        mask = bids(root='results/diffparc',template='{template}',label='{seed}',suffix='mask.nii.gz'),
        probtrack_dir = bids(root='results/diffparc',subject='{subject}',label='{seed}',space='{template}',suffix='probtrack')
    params:
        connmap_3d = lambda wildcards, input: expand(bids(root=input.probtrack_dir,include_subject_dir=False,subject=wildcards.subject,prefix='seeds_to',label='{target}',suffix='mask.nii.gz'), target=targets),
    output:
        connmap_npz = bids(root='results/diffparc',subject='{subject}',label='{seed}',space='{template}',suffix='connMap.npz')
    log: 'logs/save_connmap_to_template_npz/sub-{subject}_{seed}_{template}.log'
    group: 'participant1'
    conda: '../envs/sklearn.yml'
    script: '../scripts/save_connmap_template_npz.py'

#space-{template}
rule gather_connmap_group:
    input:
        connmap_npz = expand(bids(root='results/diffparc',subject='{subject}',label='{seed}',space='{template}',suffix='connMap.npz'), subject=subjects,allow_missing=True)
    output:
        connmap_group_npz = bids(root='results/diffparc',template='{template}',desc='concat',label='{seed}',from_='group',suffix='connMap.npz')
    log: 'logs/gather_connmap_group/{seed}_{template}.log'
    conda: '../envs/sklearn.yml'
    group: 'group1'
    script: '../scripts/gather_connmap_group.py'
     

#space-{template},  dseg
rule spectral_clustering:
    input:
        connmap_group_npz = bids(root='results/diffparc',template='{template}',desc='concat',label='{seed}',from_='group',suffix='connMap.npz')
    params:
        max_k = config['max_k']
    output:
        cluster_k = expand(bids(root='results/diffparc',template='{template}',label='{seed}',from_='group',method='spectralcosine',k='{k}',suffix='dseg.nii.gz'),k=range(2,config['max_k']+1),allow_missing=True)
    log: 'logs/spectral_clustering/{seed}_{template}.log'
    conda: '../envs/sklearn.yml'
    group: 'group1'
    script: '../scripts/spectral_clustering.py'
        
 
  

