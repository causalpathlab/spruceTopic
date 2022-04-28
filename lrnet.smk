cell_type = 'pbmc'
configfile: '/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/config/'+cell_type+'.yaml'


home_dir = config['home']
p_scripts = os.path.join(home_dir, 'sprucetopic/')
r_scripts = os.path.join(home_dir, 'scripts/')
model_file = home_dir + config['output'] + config['lr_model']['out'] + config['lr_model']['mfile']
print(home_dir)
rule all:
    input:
        expand(model_file + 'loss_plot.pdf'),
        
        expand(model_file + 'top_5_genes_topic.tsv'),
        expand(model_file + 'top_5_lrpair_topic.tsv'),
        expand(model_file + 'model_summary.csv'),

        expand(model_file + 'top5_genes_hmap.pdf'),

        expand(model_file + 'top5_lr_pair_topic_hmap.pdf'),
        expand(model_file + 'beta1_bias.pdf'),
        expand(model_file + 'beta2_bias.pdf'),

        expand(model_file + 'summary_plot.pdf')
        
rule plt_loss:
    input:
        script = r_scripts + 'fig_2_lr_loss.R'
    output:
        out_loss = model_file + 'loss_plot.pdf'
    params:
        ct = cell_type
    shell: 
        'Rscript {input.script} {params.ct}'

rule eval_model:
    input:
        script = p_scripts + 'runner.py'
    output:
        out_topic_topgenes = model_file + 'top_5_genes_topic.tsv',
        out_topic_lrpair = model_file + 'top_5_lrpair_topic.tsv',
        out_msum = model_file + 'model_summary.csv'
    params:
        mode = 'lr_net'
    shell: 
        'python {input.script} {params.mode}'

rule plt_wthmap:
    input:
        script = r_scripts + 'fig_2_lr_hmap.R',
        top_genes = rules.eval_model.output.out_topic_topgenes

    output:
        out_loss = model_file + 'top5_genes_hmap.pdf'
    params:
        ct = cell_type
    shell: 
        'Rscript {input.script} {params.ct}'

rule plt_wthmap_lrpair:
    input:
        script = r_scripts + 'fig_2_lr_hmap_tpwise.R',
        top_genes = rules.eval_model.output.out_topic_lrpair
    output:
        out_loss = model_file + 'top5_lr_pair_topic_hmap.pdf'
    params:
        ct = cell_type
    shell: 
        'Rscript {input.script} {params.ct}'

rule plt_wthmap_bias:
    input:
        script = r_scripts + 'fig_2_lr_bias_hmap.R',
        top_genes = rules.eval_model.output.out_topic_lrpair
    output:
        out_bias1 = model_file + 'beta1_bias.pdf',
        out_bias2 = model_file + 'beta2_bias.pdf'
    params:
        ct = cell_type
    shell: 
        'Rscript {input.script} {params.ct}'

rule plt_summary:
    input:
        script = r_scripts + 'fig_3_summary.R',
        top_genes = rules.eval_model.output.out_topic_lrpair
    output:
        out_summary = model_file + 'summary_plot.pdf'
    params:
        ct = cell_type
    shell: 
        'Rscript {input.script} {params.ct}'
