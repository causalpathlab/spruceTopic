configfile: '/home/BCCRC.CA/ssubedi/projects/tumour_immune_interaction/config/tcell.yaml'


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
        expand(model_file + 'marker_genes_topic.tsv'),
        
        expand(model_file + 'pp_weightmat_plot_lr.pdf'),

        expand(model_file + 'topic_lr_pair_topic.pdf')
        
rule plt_loss:
    input:
        script = r_scripts + 'fig_2_lr_loss.R'
    output:
        out_loss = model_file + 'loss_plot.pdf'
    shell: 
        'Rscript {input.script}'

rule eval_model:
    input:
        script = p_scripts + 'runner.py'
    output:
        out_topic_topgenes = model_file + 'top_5_genes_topic.tsv',
        out_topic_lrpair = model_file + 'top_5_lrpair_topic.tsv',
        out_topic_markergenes = model_file + 'marker_genes_topic.tsv'
    params:
        mode = 'lr_net'
    shell: 
        'python {input.script} {params.mode}'

rule plt_wthmap:
    input:
        script = r_scripts + 'fig_2_lr_hmap.R',
        top_genes = rules.eval_model.output.out_topic_topgenes

    output:
        out_loss = model_file + 'pp_weightmat_plot_lr.pdf'
    shell: 
        'Rscript {input.script}'

rule plt_wthmap_lrpair:
    input:
        script = r_scripts + 'fig_2_lr_hmap_tpwise.R',
        top_genes = rules.eval_model.output.out_topic_lrpair

    output:
        out_loss = model_file + 'topic_lr_pair_topic.pdf'
    shell: 
        'Rscript {input.script}'
