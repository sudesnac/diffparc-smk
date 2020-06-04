
rule combine_lr_hcp:
    input:
        lh = hcp_mmp_to_native('results/hcp_mmp/sub-{subject}/lh.native.hcp-mmp.nii.gz'),
        rh = hcp_mmp_to_native('results/hcp_mmp/sub-{subject}/rh.native.hcp-mmp.nii.gz')
    output:
        lh_rh = 'results/diffparc/sub-{subject}/masks/lh_rh.native.hcp-mmp.nii.gz'
    container: config['singularity_neuroglia']
    shell: 'fslmaths {input.lh} -max {input.rh} {output.lh_rh}'


rule get_template_seed:
    input: 
        seed = lambda wildcards: config['template_prob_seg'][wildcards.seed],
    output: 'results/template_masks/sub-{template}_desc-{seed}_probseg.nii.gz'
    singularity: config['singularity_neuroglia']
    log: 'logs/get_template_seed/{template}_{seed}.log'
    shell:
        'cp -v {input} {output} &> {log}'
 


#transform probabilistic seed to subject
rule transform_to_subject:
    input: 
        seed = rules.get_template_seed.output,
        affine =  config['ants_affine_mat'],
        invwarp =  config['ants_invwarp_nii'],
        ref = 'results/diffparc/sub-{subject}/masks/lh_rh.native.hcp-mmp.nii.gz'
    output: 'results/diffparc/sub-{subject}/masks/seed_from-{template}_{seed}.nii.gz'
    envmodules: 'ants'
    singularity: config['singularity_neuroglia']
    log: 'logs/transform_to_subject/{template}_sub-{subject}_{seed}.log'
    shell:
        'antsApplyTransforms -d 3 --interpolation Linear -i {input.seed} -o {output} -r {input.ref} -t [{input.affine},1] -t {input.invwarp} &> {log}'


    
rule resample_targets:
    input: 
        dwi = join(config['fsl_bedpost_dir'],config['bedpost_mean_s0_samples']),
        targets = 'results/diffparc/sub-{subject}/masks/lh_rh.native.hcp-mmp.nii.gz'
    params:
        seed_resolution = config['probtrack']['seed_resolution']
    output:
        mask = 'results/diffparc/sub-{subject}/masks/brain_mask_dwi.nii.gz',
        mask_res = 'results/diffparc/sub-{subject}/masks/brain_mask_dwi_resampled.nii.gz',
        targets_res = 'results/diffparc/sub-{subject}/masks/lh_rh_targets_dwi.nii.gz'
    singularity: config['singularity_neuroglia']
    log: 'logs/resample_targets/sub-{subject}.log'
    group: 'pre_track'
    shell:
        'fslmaths {input.dwi} -bin {output.mask} &&'
        'mri_convert {output.mask} -vs {params.seed_resolution} {params.seed_resolution} {params.seed_resolution} {output.mask_res} -rt nearest &&'
        'reg_resample -flo {input.targets} -res {output.targets_res} -ref {output.mask_res} -NN 0  &> {log}'



rule resample_seed:
    input: 
        seed = rules.transform_to_subject.output,
        mask_res = 'results/diffparc/sub-{subject}/masks/brain_mask_dwi_resampled.nii.gz'
    output:
        seed_res = 'results/diffparc/sub-{subject}/masks/seed_from-{template}_{seed}_resampled.nii.gz',
    singularity: config['singularity_neuroglia']
    log: 'logs/resample_seed/{template}_sub-{subject}_{seed}.log'
    group: 'pre_track'
    shell:
        #linear interp here now, since probabilistic seg
        'reg_resample -flo {input.seed} -res {output.seed_res} -ref {input.mask_res}  &> {log}'

rule binarize_subject_seed:
    input: rules.resample_seed.output
    params:
        threshold = config['prob_seg_threshold']
    output: 'results/diffparc/sub-{subject}/masks/seed_from-{template}_{seed}.binary.nii.gz'
    singularity: config['singularity_neuroglia']
    log: 'logs/binarize_subject_seed/{template}_sub-{subject}_{seed}.log'
    container: config['singularity_neuroglia']
    shell:
        'fslmaths {input} -thr {params.threshold} {output} &> {log}'
        


 

rule split_targets:
    input: 
        targets = 'results/diffparc/sub-{subject}/masks/lh_rh_targets_dwi.nii.gz',
    params:
        target_nums = lambda wildcards: [str(i) for i in range(len(targets))],
        target_seg = expand('results/diffparc/sub-{subject}/targets/{target}.nii.gz',target=targets,allow_missing=True)
    output:
        target_seg_dir = directory('results/diffparc/sub-{subject}/targets')
    singularity: config['singularity_neuroglia']
    log: 'logs/split_targets/sub-{subject}.log'
    threads: 32 
    group: 'pre_track'
    shell:
        'mkdir -p {output} && parallel  --jobs {threads} fslmaths {input.targets} -thr {{1}} -uthr {{1}} -bin {{2}} &> {log} ::: {params.target_nums} :::+ {params.target_seg}'

rule gen_targets_txt:
    input:
        target_seg_dir = 'results/diffparc/sub-{subject}/targets'
    params:
        target_seg = expand('results/diffparc/sub-{subject}/targets/{target}.nii.gz',target=targets,allow_missing=True)
    output:
        target_txt = 'results/diffparc/sub-{subject}/target_images.txt'
    log: 'logs/get_targets_txt/sub-{subject}.log'
    group: 'pre_track'
    run:
        f = open(output.target_txt,'w')
        for s in params.target_seg:
            f.write(f'{s}\n')
        f.close()


rule run_probtrack:
    input:
        seed_res = rules.binarize_subject_seed.output,
        target_txt = rules.gen_targets_txt.output,
        mask = 'results/diffparc/sub-{subject}/masks/brain_mask_dwi.nii.gz',
        target_seg_dir = 'results/diffparc/sub-{subject}/targets'
    params:
        bedpost_merged = join(config['fsl_bedpost_dir'],config['bedpost_merged_prefix']),
        probtrack_opts = config['probtrack']['opts'],
        out_target_seg = expand('results/diffparc/sub-{subject}/probtrack_{template}_{seed}/seeds_to_{target}.nii.gz',target=targets,allow_missing=True)
    output:
        probtrack_dir = directory('results/diffparc/sub-{subject}/probtrack_{template}_{seed}')
    threads: 2
    resources: 
        mem_mb = 8000, 
        time = 30, #30 mins
        gpus = 1 #1 gpu
    log: 'logs/run_probtrack/{template}_sub-{subject}_{seed}.log'
    shell:
        'mkdir -p {output.probtrack_dir} && probtrackx2_gpu --samples={params.bedpost_merged}  --mask={input.mask} --seed={input.seed_res} ' 
        '--targetmasks={input.target_txt} --seedref={input.seed_res} --nsamples={config[''probtrack''][''nsamples'']} ' 
        '--dir={output.probtrack_dir} {params.probtrack_opts} -V 2  &> {log}'


rule transform_conn_to_template:
    input:
        connmap_dir = 'results/diffparc/sub-{subject}/probtrack_{template}_{seed}',
        affine =  config['ants_affine_mat'],
        warp =  config['ants_warp_nii'],
        ref = config['ants_ref_nii']
    params:
        in_connmap_3d = expand('results/diffparc/sub-{subject}/probtrack_{template}_{seed}/seeds_to_{target}.nii.gz',target=targets,allow_missing=True),
        out_connmap_3d = expand('results/diffparc/sub-{subject}/probtrack_{template}_{seed}_warped/seeds_to_{target}_space-{template}.nii.gz',target=targets,allow_missing=True)
    output:
        connmap_dir = directory('results/diffparc/sub-{subject}/probtrack_{template}_{seed}_warped')
    envmodules: 'ants'
    singularity: config['singularity_neuroglia']
    threads: 32
    resources:
        mem_mb = 128000
    log: 'logs/transform_conn_to_template/sub-{subject}_{seed}_{template}.log'
    group: 'post_track'
    shell:
        'mkdir -p {output} && ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1 parallel  --jobs {threads} antsApplyTransforms -d 3 --interpolation Linear -i {{1}} -o {{2}}  -r {input.ref} -t {input.warp} -t {input.affine} &> {log} :::  {params.in_connmap_3d} :::+ {params.out_connmap_3d}' 



rule binarize_template_seed:
    input: rules.get_template_seed.output
    params:
        threshold = config['prob_seg_threshold']
    output: 'results/template_masks/sub-{template}_desc-{seed}_mask.nii.gz'
    singularity: config['singularity_neuroglia']
    log: 'logs/binarize_template_seed/{template}_{seed}.log'
    shell:
        'fslmaths {input} -thr {params.threshold} {output} &> {log}'
        



rule save_connmap_template_npz:
    input:
        mask = 'results/template_masks/sub-{template}_desc-{seed}_mask.nii.gz',
        connmap_dir = 'results/diffparc/sub-{subject}/probtrack_{template}_{seed}_warped'
    params:
        connmap_3d = expand('results/diffparc/sub-{subject}/probtrack_{template}_{seed}_warped/seeds_to_{target}_space-{template}.nii.gz',target=targets,allow_missing=True),
    output:
        connmap_npz = 'results/diffparc/sub-{subject}/connmap/sub-{subject}_space-{template}_seed-{seed}_connMap.npz'
    log: 'logs/save_connmap_to_template_npz/sub-{subject}_{seed}_{template}.log'
    group: 'post_track'
    script: '../scripts/save_connmap_template_npz.py'

rule gather_connmap_group:
    input:
        connmap_npz = expand('results/diffparc/sub-{subject}/connmap/sub-{subject}_space-{template}_seed-{seed}_connMap.npz',subject=subjects,allow_missing=True)
    output:
        connmap_group_npz = 'results/connmap/group_space-{template}_seed-{seed}_connMap.npz'
    log: 'logs/gather_connmap_group/{seed}_{template}.log'
    script: '../scripts/gather_connmap_group.py'
     
rule spectral_clustering:
    input:
        connmap_group_npz = 'results/connmap/group_space-{template}_seed-{seed}_connMap.npz'
    params:
        max_k = config['max_k']
    output:
        cluster_k = expand('results/clustering/group_space-{template}_seed-{seed}_method-spectralcosine_k-{k}_cluslabels.nii.gz',k=range(2,config['max_k']+1),allow_missing=True)
    script: '../scripts/spectral_clustering.py'
        
     
    

    


