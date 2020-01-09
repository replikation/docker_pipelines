process sourclusterPlot {
    publishDir "${params.output}/${name}", mode: 'copy', pattern: "cluster_results.html"
    label 'bokeh'
    echo true
  input:
    tuple val(name), file(results) 
  output:
	  tuple val(name), file("cluster_results.html")
  script:
    """
    #!/usr/bin/env python

    from collections import defaultdict

    import numpy as np
    import scipy.cluster.hierarchy as sch

    from bokeh.plotting import figure, show, save, output_file
    from bokeh.sampledata.les_mis import data

    with open("${results}", 'r') as file: 
        header = next(file).strip().split(',')
        D = []  # distance matrix
        for source in header:
            values = [float(x) for x in next(file).strip().split(',')]
            D.append(values)

    D_ = np.array(D)
    Y = sch.linkage(D_, method='single')
    Z = sch.dendrogram(Y, orientation='right')
    ix = Z['leaves']

    D_ = D_[ix, :]
    D_ = D_[:, ix]
    header_sorted = [header[i] for i in ix]

    l = []
    for i_, i in enumerate(D_):
        for j_, j in enumerate(i):
            l.append([header_sorted[i_], header_sorted[j_], j, '#225ea8'])   # l.append([header_sorted[i_], header_sorted[j_], j, '#ff7f00'])

    xname, yname, alpha, color = zip(*l)

    data=dict(
        xname=xname,
        yname=yname,
        colors=color,
        alphas=alpha,
        )

    p = figure(title='Clustered Jaccard distance (sourmash)',
              x_axis_location='above', tools='hover,save',
              x_range=list(reversed(header_sorted)), y_range=header_sorted,
              tooltips = [('names', '@yname, @xname'), ('distance', '@alphas')]
              )

    p.plot_width = ${params.size}
    p.plot_height = ${params.size}
    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = '5pt'
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = np.pi/3

    p.rect('xname', 'yname', 0.9, 0.9, source=data,
          color='colors', alpha='alphas', line_color=None,
          hover_line_color='black', hover_color='colors')


    output_file("cluster_results.html", title='Clustered Jaccard distance (sourmash)')
    save(p)
    """
}

// tempDir may be neccessary
/*

    code by @phiweger

    A :math:`(n-1)` by 4 matrix ``Z`` is returned. At the
    :math:`i`-th iteration, clusters with indices ``Z[i, 0]`` and
    ``Z[i, 1]`` are combined to form cluster :math:`n + i`. A
    cluster with an index less than :math:`n` corresponds to one of
    the :math:`n` original observations. The distance between
    clusters ``Z[i, 0]`` and ``Z[i, 1]`` is given by ``Z[i, 2]``. The
    fourth value ``Z[i, 3]`` represents the number of original
    observations in the newly formed cluster.

    # How to sort a distance matrix
    # https://gmarti.gitlab.io/ml/2017/09/07/how-to-sort-distance-matrix.html
    # THIS led me to the right solution -- stackoverflow.com/questions/2455761/
    # Reorder, first rows, then columns

    # https://github.com/dib-lab/sourmash/blob/master/sourmash/commands.py#L235

    
    # colormap = ['#444444', '#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99',
    #             '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a']

        #parser = argparse.ArgumentParser(description='Process some integers.')
    #parser.add_argument('-d', help='Distance matrix from sourmash')
    #parser.add_argument('-o', help='File to save the plot to')
    #args = parser.parse_args()
*/