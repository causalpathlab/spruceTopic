configfile: '/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/config/pbmc.yaml'


home_dir = config['home']
p_scripts = os.path.join(home_dir, 'sprucetopic/')
r_scripts = os.path.join(home_dir, 'scripts/')
model_file = home_dir + config['output'] + config['nbr_model']['out'] + config['nbr_model']['mfile']

rule all:
    input:
        expand(model_file + 'loss_plot.pdf'),
        
        expand(model_file + 'umap_zz_topic_label.png'),
        expand(model_file + 'umap_zz_scanpy_cluster_label.png'),
        expand(model_file + 'hh_cell_topic_sample.tsv'),
        expand(model_file + 'top_5_genes_topic.tsv'),

        expand(model_file + 'cell_topic_struct_plot_pbmc.pdf'),
        
        expand(model_file + 'topic_gene_weight_hmap.pdf'),

        expand(model_file + 'marker_gene_weight_hmap.pdf')
        


rule plt_loss:
    input:
        script = r_scripts + 'fig_1_loss.R'
    output:
        out_loss = model_file + 'loss_plot.pdf'
    shell: 
        'Rscript {input.script}'

rule eval_model:
    input:
        script = p_scripts + 'runner.py'
    output:
        out_umap_topic = model_file + 'umap_zz_topic_label.png',
        out_umap_topic_celllabel = model_file + 'umap_zz_scanpy_cluster_label.png',
        out_sample_topic = model_file + 'hh_cell_topic_sample.tsv',
        out_topic_topgenes = model_file + 'top_5_genes_topic.tsv',
        out_topic_markergenes = model_file + 'marker_genes_topic.tsv'
    params:
        mode = 'nbr_net'
    shell: 
        'python {input.script} {params.mode}'

rule plt_st:
    input:
        script = r_scripts + 'fig_1_stplot.R',
        latent = rules.eval_model.output.out_sample_topic
    output:
        out_loss = model_file + 'cell_topic_struct_plot_pbmc.pdf'
    shell: 
        'Rscript {input.script}'

rule plt_wthmap:
    input:
        script = r_scripts + 'fig_1_hmap.R',
        top_genes = rules.eval_model.output.out_topic_topgenes

    output:
        out_loss = model_file + 'topic_gene_weight_hmap.pdf'
    shell: 
        'Rscript {input.script}'

rule plt_markerg_hmap:
    input:
        script = r_scripts + 'fig_1_hmap_mg.R',
        top_mgenes = rules.eval_model.output.out_topic_markergenes

    output:
        out_loss = model_file + 'marker_gene_weight_hmap.pdf'
    shell: 
        'Rscript {input.script}'
