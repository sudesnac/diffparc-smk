#warp group-based clusters back to each subject
#space-T1w, desc-groupclus?, dseg
rule transform_clus_to_subj:
    input: 
        cluster_k = expand(bids(root='results/diffparc',template='{template}',label='{seed}',from_='group',method='spectralcosine',k='{k}',suffix='dseg.nii.gz'),k=range(2,config['max_k']+1),allow_missing=True),
        affine =  config['ants_affine_mat'],
        invwarp =  config['ants_invwarp_nii'],
        ref = bids(root='results/diffparc',subject='{subject}',space='individual',label=config['targets_atlas_name'],suffix='dseg.nii.gz')
    output: 
        cluster_k = expand(bids(root='results/tractmap',subject='{subject}',space='individual',label='{seed}',method='spectralcosine',k='{k}',from_='{template}',suffix='dseg.nii.gz'),k=range(2,config['max_k']+1),allow_missing=True)
    envmodules: 'ants'
    container: config['singularity']['neuroglia']
    log: 'logs/transform_clus_to_subject/sub-{subject}_template-{template}_{seed}.log'
    group: 'participant2'
    threads: 8
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1 parallel  --jobs {threads} antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {{1}} -o {{2}}  -r {input.ref} -t [{input.affine},1] -t {input.invwarp} &> {log} :::  {input.cluster_k} :::+ {output.cluster_k}' 



#create brainmask from bedpost data, and resample to chosen resolution
#space-T1w res-? mask
rule resample_brainmask_tractmaps:
    input: 
        dwi = join(config['fsl_bedpost_dir'],config['bedpost_mean_s0_samples']),
    params:
        seed_resolution = config['probtrack_tractmap']['seed_resolution']
    output:
        mask = bids(root='results/tractmap',subject='{subject}',label='brain',suffix='mask.nii.gz'),
        mask_res = bids(root='results/tractmap',subject='{subject}',label='brain',res='super',suffix='mask.nii.gz'),
    container: config['singularity']['neuroglia']
    log: 'logs/resample_brainmask_tractmaps/sub-{subject}.log'
    group: 'participant2'
    shell:
        'fslmaths {input.dwi} -bin {output.mask} &&'
        'mri_convert {output.mask} -vs {params.seed_resolution} {params.seed_resolution} {params.seed_resolution} {output.mask_res} -rt nearest &> {log}'


#resamples clus seg to dwi brainmask resolution
#space-T1w, res=?   dseg
rule resample_clus_seed:
    input: 
        seed = bids(root='results/tractmap',subject='{subject}',space='individual',label='{seed}',method='spectralcosine',k='{k}',from_='{template}',suffix='dseg.nii.gz'),
        mask_res = bids(root='results/tractmap',subject='{subject}',label='brain',res='super',suffix='mask.nii.gz'),
    output:
        seed_res = bids(root='results/tractmap',subject='{subject}',space='individual',label='{seed}',method='spectralcosine',k='{k}',from_='{template}',res='super',suffix='dseg.nii.gz'),
    container: config['singularity']['neuroglia']
    log: 'logs/resample_clus_seed/sub-{subject}_template-{template}_{seed}_k-{k}.log'
    group: 'participant2'
    shell:
        #linear interp here now, since probabilistic seg
        'reg_resample -flo {input.seed} -res {output.seed_res} -ref {input.mask_res} -NN 0 &> {log}'



#split segmentation into binary masks
# space-T1w   mask
rule subj_split_clus_to_binary_masks:
    input: 
        cluster_k = bids(root='results/tractmap',subject='{subject}',space='individual',label='{seed}',method='spectralcosine',k='{k}',from_='{template}',res='super',suffix='dseg.nii.gz'),
    params:
        mask_file = lambda wildcards, output: bids(root=output.cluster_k_splitdir,subject=wildcards.subject,label='%02d',suffix='mask.nii.gz',include_subject_dir=False),
        mask_bg = lambda wildcards, output: bids(root=output.cluster_k_splitdir,subject=wildcards.subject,label='00',suffix='mask.nii.gz',include_subject_dir=False) #we remove this file 
    output:
        cluster_k_splitdir = directory(bids(root='results/tractmap',subject='{subject}',space='individual',label='{seed}',method='spectralcosine',k='{k}',from_='{template}',res='super',suffix='seeds'))
    container: config['singularity']['neuroglia']
    log: 'logs/subj_split_clus_to_binary_masks/sub-{subject}_template-{template}_{seed}_k-{k}.log'
    group: 'participant2'
    threads: 8
    shell:
        #use c3d's split command to go from discrete seg image to multiple binary images.. we remove image 00 since that is the background label
        'mkdir {output.cluster_k_splitdir} && c3d {input.cluster_k} -split -oo {params.mask_file} &>> {log}  && rm -f {params.mask_bg}'


#perform tracking from each cluster in subj space to get tract maps
# check bids-deriv dwi - tractography type?
# space-T1w, res-?
rule track_from_clusters:
    input:
        cluster_k_splitdir = bids(root='results/tractmap',subject='{subject}',space='individual',label='{seed}',method='spectralcosine',k='{k}',from_='{template}',res='super',suffix='seeds'),
        mask = bids(root='results/tractmap',subject='{subject}',label='brain',suffix='mask.nii.gz'),
    params:
        seeds = lambda wildcards, input: expand(bids(root=input.cluster_k_splitdir,subject=wildcards.subject,label='{k_index:02d}',suffix='mask.nii.gz',include_subject_dir=False), k_index=range(1,int(wildcards.k)+1)),
        bedpost_merged = join(config['fsl_bedpost_dir'],config['bedpost_merged_prefix']),
        probtrack_opts = config['probtrack_tractmap']['opts'],
        nsamples = config['probtrack_tractmap']['nsamples'],
        out_track_dirs = lambda wildcards,output: expand('{}/label-{{k_index:02d}}'.format(output.probtrack_dir),k_index=range(1,int(wildcards.k)+1)),
        container = '/project/6050199/akhanf/singularity/bids-apps/fsl_6.0.3_cuda9.1.sif'
    output:
        probtrack_dir = directory(bids(root='results/tractmap',subject='{subject}',space='individual',label='{seed}',method='spectralcosine',k='{k}',from_='{template}',res='super',suffix='probtrack')),
    threads: 8
    resources:
        mem_mb = 32000,
        time = 30, #30 mins
        gpus = 1 #1 gpu
    log: 'logs/track_from_clusters/sub-{subject}_template-{template}_{seed}_k-{k}.log'
    group: 'participant2'
    shell:
        #this job runs probtrack for each seed
        'mkdir -p {params.out_track_dirs} && '
        ' parallel --jobs 1 '
        '   singularity exec -e --nv {params.container} probtrackx2_gpu --samples={params.bedpost_merged}  --mask={input.mask} --seed={{1}} '
        '    --seedref={{1}} --nsamples={params.nsamples} '
        '    --dir={{2}} {params.probtrack_opts} -V 2  &> {log}  '
        ' ::: {params.seeds} :::+ {params.out_track_dirs}'

# check bids-deriv dwi - tractography ?
# space-T1w, res-?
rule combine_tractmaps:
    input:
        probtrack_dir = bids(root='results/tractmap',subject='{subject}',space='individual',label='{seed}',method='spectralcosine',k='{k}',from_='{template}',res='super',suffix='probtrack')
    params:
        tractmaps = lambda wildcards,input: expand('{}/label-{{k_index:02d}}/fdt_paths.nii.gz'.format(input.probtrack_dir),k_index=range(1,int(wildcards.k)+1)),
    output:
        tractmaps_4d = bids(root='results/tractmap',subject='{subject}',space='individual',label='{seed}',method='spectralcosine',k='{k}',from_='{template}',res='super',suffix='tractmap4d.nii.gz')
    log: 'logs/combine_tractmaps/sub-{subject}_template-{template}_{seed}_k-{k}.log'
    container: config['singularity']['neuroglia']
    group: 'participant2'
    resources:
        mem_mb = 32000 
    shell:
        'fslmerge -t {output.tractmaps_4d} {params.tractmaps} &> {log}'


#transform tract maps back to template space
# space-{template}, tractogrpahy?
rule transform_tractmaps_to_template:
    input:
        probtrack_dir = bids(root='results/tractmap',subject='{subject}',space='individual',label='{seed}',method='spectralcosine',k='{k}',from_='{template}',res='super',suffix='probtrack'),
        affine =  config['ants_affine_mat'],
        warp =  config['ants_warp_nii'],
        ref = config['ants_ref_nii']
    params:
        in_tractmaps = lambda wildcards,input: expand('{}/label-{{k_index:02d}}/fdt_paths.nii.gz'.format(input.probtrack_dir),k_index=range(1,int(wildcards.k)+1)),
        out_tractmaps = lambda wildcards,output: expand('{}/label-{{k_index:02d}}_tractmap.nii.gz'.format(output.probtrack_dir),k_index=range(1,int(wildcards.k)+1))
    output:
        probtrack_dir = directory(bids(root='results/tractmap',subject='{subject}',label='{seed}',method='spectralcosine',k='{k}',space='{template}',suffix='probtrack')),
    container: config['singularity']['neuroglia']
    log: 'logs/transform_tractmaps_to_template/sub-{subject}_{seed}_{template}_k-{k}.log'
    group: 'participant2'
    threads: 8
    resources:
        mem_mb = 32000,
        time = 120
    shell:
        'mkdir -p {output} && ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1 parallel  --jobs {threads} antsApplyTransforms -d 3 --interpolation Linear -i {{1}} -o {{2}}  -r {input.ref} -t {input.warp} -t {input.affine} &> {log} :::  {params.in_tractmaps} :::+ {params.out_tractmaps}' 


# space-{template}, tractography 4d?
rule combine_tractmaps_warped:
    input:
        probtrack_dir = bids(root='results/tractmap',subject='{subject}',label='{seed}',method='spectralcosine',k='{k}',space='{template}',suffix='probtrack')
    params:
        tractmaps = lambda wildcards,input: expand('{}/label-{{k_index:02d}}_tractmap.nii.gz'.format(input.probtrack_dir),k_index=range(1,int(wildcards.k)+1))
    output:
        tractmaps_4d = bids(root='results/tractmap',subject='{subject}',label='{seed}',method='spectralcosine',k='{k}',space='{template}',suffix='tractmap4d.nii.gz')
    log: 'logs/combine_tractmaps_warped/sub-{subject}_template-{template}_{seed}_k-{k}.log'
    container: config['singularity']['neuroglia']
    resources:
        mem_mb = 32000 
    group: 'participant2'
    shell:
        'fslmerge -t {output.tractmaps_4d} {params.tractmaps} &> {log}'



#now, average the tract maps (over subjects) in template space
# space-template, desc-avg , tractography 4d?
rule avg_tractmaps_template:
    input: 
        tractmaps_4d = expand(bids(root='results/tractmap',subject='{subject}',label='{seed}',method='spectralcosine',k='{k}',space='{template}',suffix='tractmap4d.nii.gz'),subject=subjects,allow_missing=True),
    output: 
        average = bids(root='results/tractmap',template='{template}',label='{seed}',method='spectralcosine',k='{k}',desc='average',suffix='tractmap4d.nii.gz')
    container: config['singularity']['neuroglia']
    group: 'group2'
    threads: 8
    resources:
        mem_mb = 32000 
    shell:
        'AverageImages 4 {output} 0 {input}'
             

#use voting to get a discrete segmentation of tract map 
# space-template, desc-avgtractmap, dseg
rule vote_tractmap_template:
    input: 
        tractmaps = bids(root='results/tractmap',template='{template}',label='{seed}',method='spectralcosine',k='{k}',desc='average',suffix='tractmap4d.nii.gz')
    params:
        bg_th = 500   # set only if avg streamline count > bg_th 
    output:
        vote_seg = bids(root='results/tractmap',template='{template}',label='{seed}',method='spectralcosine',k='{k}',desc='avgtractmap',suffix='dseg.nii.gz')
    container: config['singularity']['neuroglia']
    log: 'logs/vote_tractmap_template/{template}_{seed}_{k}.log'
    group: 'group2'
    threads: 8
    resources:
        mem_mb = 32000
    shell: 
        #get first vol, set it to bg_th value, this becomes image 0
        # then load all the tractmaps as images 1-k
        # then vote amongst the 0-k images, those where bg (bg_th) is highest get set to 0, otherwise set to the index (from 1-k)
        'c4d -verbose {input.tractmaps} -slice w 0:0 -threshold -inf +inf {params.bg_th} 0 {input.tractmaps} -slice w 0:-1  -vote -type uchar {output.vote_seg} &> {log}' 


