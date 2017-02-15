#!/usr/bin/env python

'''NB! matplotlib 1.5.2 doesn't work with 80x workspace! see more at https://github.com/cms-sw/cmssw/issues/15660
       Therefore we must tweak the PYTHONPATH variable a little bit so that we'd use matplotlib 1.2.1
'''

import sys, os
sys.path = ['/cvmfs/cms.cern.ch/%s/external/py2-matplotlib/1.2.1-ikhhed/lib/python2.7/site-packages/' % os.getenv('SCRAM_ARCH')] + sys.path

import argparse, logging, matplotlib.pyplot as plt, matplotlib.ticker as tck, numpy as np

def auc(x, y):
  '''Calculates area under curve (AUC)
  :param x: float array, x-coordinates of the curve
  :param y: float array, y-coordinates of the curve
  :return: float, AUC
  '''
  s = 0.
  for i in range(len(x) - 1):
    s += (x[i] - x[i + 1]) * (y[i + 1] + y[i]) / 2.
  return s

if __name__ == '__main__':
  logging.basicConfig(
    stream = sys.stdout,
    level  = logging.INFO,
    format = '%(asctime)s - %(levelname)s: %(message)s'
  )

  # set the help description width to 45 characters
  class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
      if text.startswith('R|'):
        return text[2:].splitlines()
      return argparse.HelpFormatter._split_lines(self, text, width)

  parser = argparse.ArgumentParser(formatter_class = lambda prog: SmartFormatter(prog, max_help_position = 50))
  parser.add_argument('-i', '--input', metavar = 'directory', required = True, type = str, default = '',
                      help = 'R|Base directory where the MEM results are stored')
  parser.add_argument('-l', '--list', metavar = 'dir', required = False, type = str, default = 'mem', nargs = '+',
                      help = 'R|List of subdirectories containing MEM results')
  parser.add_argument('-t', '--title', metavar = 'title', required = False, type = str, default = 'MEM ROC curve',
                      help = 'R|Title of the plot')
  parser.add_argument('-x', '--xlabel', metavar = 'text', required = False, type = str, default = 'ttHJetToNonbb_M125',
                      help = 'R|Label on the x-axis')
  parser.add_argument('-y', '--ylabel', metavar = 'text', required = False, type = str, default = 'TTZToLLNuNu',
                      help = 'R|Label on the y-axis')
  parser.add_argument('-b', '--branch-name', metavar = 'text', required = False, type = str, default = 'lhRatioNP',
                      help = 'R|Name of the likelihood ratio branch')
  parser.add_argument('-o', '--output', metavar = 'file', required = False, type = str, default = '',
                      help = 'R|Output file')
  parser.add_argument('-e', '--error-bands', dest = 'error_bands', action = 'store_true', default = False,
                      help = 'R|Plot error bands around the ROC curves')
  parser.add_argument('-f', '--force', dest = 'force', action = 'store_true', default = False,
                      help = 'R|Force the creation of output directory if missing')
  parser.add_argument('-v', '--verbose', dest = 'verbose', action = 'store_true', default = False,
                      help = 'R|Enable verbose printout')
  args = parser.parse_args()

  if args.verbose:
    logging.getLogger().setLevel(logging.DEBUG)

  if not os.path.isdir(args.input):
    logging.error("No such directory: {dirname}".format(dirname = args.input))
    sys.exit(1)

  if not args.output:
    plot_dir = os.path.join(args.input, 'plots')
    outputfile = os.path.join(plot_dir, '%s_%s.pdf' % (args.xlabel, args.ylabel))
  else:
    outputfile = args.output

  # testing if the list of MEM subdirs are really there
  mem_subdirs = { x : { 'dir' : os.path.join(args.input, x) } for x in args.list }
  mem_subdirs_missing = [ v['dir'] for v in mem_subdirs.values() if not os.path.exists(v['dir']) ]
  if mem_subdirs_missing:
    logging.error("The following subdirectories are missing: {subdir_list}".format(
      subdir_list = ', '.join(mem_subdirs_missing)
    ))
    sys.exit(1)

  # let's read the comment files (if they're there); consider only the first line
  # if the comment file is missing then we just use the subfolder name as the comment
  for key, mem_subdir in mem_subdirs.iteritems():
    logging.debug("Reading comment files and reassuring the existence of ROC CSV files")
    comment_file = os.path.join(mem_subdir['dir'], 'comment.txt')
    if os.path.exists(comment_file):
      comment_lines = open(comment_file, 'r').readlines()
      if comment_lines:
        mem_subdir['comment'] = comment_lines[0]
      else:
        mem_subdir['comment'] = key
    else:
      mem_subdir['comment'] = key

    roc_file = os.path.join(
      mem_subdir['dir'], 'roc', 'csvs', 'roc_{xlabel}_{ylabel}_{branch_name}.csv'.format(
        xlabel      = args.xlabel,
        ylabel      = args.ylabel,
        branch_name = args.branch_name,
      )
    )
    if not os.path.isfile(roc_file):
      logging.error("File '{roc_filename}' does not exist (check if the x- and y-labels are correct".format(
        roc_filename = roc_file,
      ))
      sys.exit(1)
    mem_subdir['roc_csv'] = roc_file

  # check if the directory to output file exists
  # we've postponed this check at latest possible b/c other checks might fail and we don't want to
  # create an empty directory if nothing won't be written there
  outputfile_parentdir = os.path.dirname(outputfile)
  if not os.path.isdir(outputfile_parentdir):
    if args.force:
      logging.debug("Directory '{output_dirname}' of the output file '{output_filename}' does not exist; "
                    "attempting to create one".format(
        output_dirname  = outputfile_parentdir,
        output_filename = outputfile,
      ))
      try:
        os.makedirs(outputfile_parentdir)
      except IOError as err:
        logging.error("Could not create directory '{output_dirname}' because: {reasons}".format(
          output_dirname = outputfile_parentdir,
          reasons        = err,
        ))
    else:
      logging.error("Directory '{output_dirname}' of the output file '{output_filename}' does not exist".format(
        output_dirname  = outputfile_parentdir,
        output_filename = outputfile,
      ))
      sys.exit(1)

  # now we're ready to plot the ROC curve
  logging.debug("Plotting ...")
  fig, ax = plt.subplots(figsize = (12, 8))
  for _, values in mem_subdirs.iteritems():
    data = np.loadtxt(values['roc_csv'], delimiter = ',', unpack = True)
    p = plt.plot(data[0], data[1], label = '%s (%.2f)' % (values['comment'], auc(data[0], data[1])), lw = 2)
    minIdx = np.argmin(np.sqrt((data[0] - 1) ** 2 + data[1] ** 2))
    xmin = data[0][minIdx]
    ymin = data[1][minIdx]
    plt.plot((xmin,), (ymin,), marker = 'o', markersize = 10, ls = '', c = 'none')
    if not np.allclose(data[3], data[5]):
      ax.fill_between(
        data[0], data[3], data[5],
        interpolate = True,
        alpha       = 0.1,
        linewidth   = 0.0,
        color       = p[0].get_color(),
      )
    plt.annotate('%.2f' % data[2][minIdx], xy = (xmin + 0.025, ymin))
    plt.annotate('Optimal cutoff points are circled with\ncorresponding WP values next to the circle',
                 xy = (0.01, 0.6), fontsize = 12, family = 'sans-serif')
  plt.plot((0., 1.), (0., 1.), '--', c = 'black', label = 'random (0.50)')
  plt.xticks(np.arange(0., 1.1, 0.1))
  plt.yticks(np.arange(0., 1.1, 0.1))
  plt.legend(loc = 2)
  plt.grid(True)
  plt.tick_params(labelright = True, labeltop = True)
  ax.xaxis.set_major_locator(tck.MultipleLocator(0.1))
  ax.xaxis.set_minor_locator(tck.MultipleLocator(0.01))
  ax.yaxis.set_major_locator(tck.MultipleLocator(0.1))
  ax.yaxis.set_minor_locator(tck.MultipleLocator(0.01))
  plt.xlabel('%s efficiency' % args.xlabel, fontsize = 14)
  plt.ylabel('%s efficiency' % args.ylabel, fontsize = 14)
  plt.suptitle(args.title, fontsize = 16)
  plt.savefig(outputfile, bbox_inches = 'tight')
  logging.info("Saved the figure to: {image_path}".format(image_path = outputfile))
