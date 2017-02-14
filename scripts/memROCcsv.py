#!/usr/bin/env python

import argparse, logging, sys, os, ROOT, array, numpy as np, itertools

def get_roc(wps, lrs):
  counts = []
  for wp in wps:
    counts.append(sum((lrs >= wp) * 1.))
  counts = np.asarray(counts)
  counts /= len(lrs)
  if 0. not in counts:
    counts[-1] = 0.
  return counts

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

  parser = argparse.ArgumentParser(formatter_class = lambda prog: SmartFormatter(prog, max_help_position = 45))
  parser.add_argument('-i', '--input', metavar = 'files', required = True, type = str, nargs = '+', default = [],
                      help = 'R|Comma-separated list of input files')
  parser.add_argument('-c', '--class-labels', metavar = 'label', required = True, type = str, nargs = '+', default = [],
                      help = 'R|Class of each input file (either \'sig\' for signal or \'bkg\' background), '
                             'as a comma-separated list')
  parser.add_argument('-s', '--signal', metavar = 'name', required = False, type = str, default = '',
                      help = 'R|Name of the signal hypothesis branch')
  parser.add_argument('-b', '--bkg', metavar = 'name', required = False, type = str, default = '',
                      help = 'R|Name of the background hypothesis branch')
  parser.add_argument('-l', '--likelihood-ratio', metavar = 'name', required = False, type = str, default = '',
                      help = 'R|Name of the likelihood ratio branch')
  parser.add_argument('-o', '--output', metavar = 'directory', required = True, type = str, default = '',
                      help = 'R|Output directory')
  parser.add_argument('-t', '--tree', metavar = 'name', required = False, type = str, default = 'tree',
                      help = 'R|TTree name (default: tree)')
  parser.add_argument('-f', '--force', dest = 'force', action = 'store_true', default = False,
                      help = 'R|Force the creation of output directory if missing')
  parser.add_argument('-v', '--verbose', dest = 'verbose', action = 'store_true', default = False,
                      help = 'R|Enable verbose printout')
  args = parser.parse_args()

  if args.verbose:
    logging.getLogger().setLevel(logging.DEBUG)

  # check if the number of input files matches to the number of class labels
  list_input_files = args.input
  list_labels      = args.class_labels
  if len(list_input_files) != len(list_labels):
    raise parser.error("Number of input files doesn't match with supplied number of class labels")

  # check if all the files actually exist
  missing_files = filter(lambda x: not os.path.isfile(x), list_input_files)
  if missing_files:
    raise parser.error("The following input files do not exist: {missing}".format(missing = ' '.join(missing_files)))

  # assert that all class labels are correct
  unrecognized_labels = filter(lambda x: x not in ['sig', 'bkg'], list_labels)
  if unrecognized_labels:
    raise parser.error("Unrecognized labels (use only 'sig' or 'bkg'): {list}".format(list = ' '.join(unrecognized_labels)))

  # assert that the list of labels contains exactly one signal label and at least one background label
  if list_labels.count('sig') != 1:
    raise parser.error("You must provide exactly one signal label ('sig')")
  if list_labels.count('bkg') < 1:
    raise parser.error("You must provide at least one background label ('bkg')")

  # check if either the signal-background pair or the likelihood ratio is supplied, but not both
  xor = lambda lhs, rhs: bool(lhs) ^ bool(rhs)
  branch_sig     = args.signal
  branch_bkg     = args.bkg
  branch_lr  = args.likelihood_ratio
  if xor(branch_sig, branch_bkg):
    raise parser.error("You must specify -s/--signal and -b/--bkg together, or not at all")
  if not xor(branch_sig and branch_bkg, branch_lr):
    raise parser.error("You must specify either -s/--signal and -b/--bkg, or -l/--likelihood-ratio, but not both or none")

  # check if the output directory exists; if not, check if force flag has been set; if not; raise an exception
  output_dir = args.output
  if not os.path.isdir(output_dir):
    if not args.force:
      raise parser.error("Directory {dirname} does not exists; use -f/--force to force-create it".format(
        dirname = output_dir,
      ))
    else:
      try:
        os.makedirs(output_dir)
      except IOError as err:
        raise parser.error("Caught an error while creating directory {dirname}: {reason}".format(
          dirname = output_dir,
          reason  = err,
        ))
      logging.debug("Created directory {dirname}".format(dirname = output_dir))

  # let's build the class label - file name pair
  input = {}
  for i in range(len(list_labels)):
    label = list_labels[i]
    if label not in input:
      input[label] = []
    input[label].append({'file' : list_input_files[i], 'lr' : [], 'lr_up' : [], 'lr_down' : [], })

  # let's find common prefix of the input files
  list_input_basenames = map(os.path.basename, list_input_files)
  common_prefix = ''.join(
    c[0] for c in itertools.takewhile(lambda x: all(x[0] == y for y in x), itertools.izip(*list_input_basenames))
  )
  common_suffix = ''.join(
    c[0] for c in itertools.takewhile(lambda x: all(x[0] == y for y in x),
                                      itertools.izip(*map(lambda x: x[::-1], list_input_basenames)))
  )[::-1]
  for label in input:
    for sig_bkg in input[label]:
      sig_bkg['unique_id'] = os.path.basename(sig_bkg['file'])[len(common_prefix):-len(common_suffix)]

  # now do the actual work here
  tree_name = args.tree
  for label in input:
    for sig_bkg_dict in input[label]:
      file_name = sig_bkg_dict['file']
      logging.debug("Opening {file_name} for reading".format(file_name = file_name))

      # check if the tree is present
      file_root = ROOT.TFile(file_name, 'read')
      if tree_name not in [x.GetName() for x in file_root.GetListOfKeys()]:
        raise ValueError("Invalid TTree name {tree_name}".format(tree_name = tree_name))

      # check if the provided branch names are correct
      tree = file_root.Get(tree_name)
      branch_names = set(x.GetName() for x in tree.GetListOfBranches())

      if branch_lr:
        branch_lr_err = '%s_err' % branch_lr
        missing_branches = {branch_lr, branch_lr_err} - branch_names
        if missing_branches:
          raise ValueError("No such branches in {tree_name}: {missing_branches}".format(
            tree_name        = tree_name,
            missing_branches = ' '.join(missing_branches)
          ))

        # read the branches, assuming type double
        lr     = array.array('d', [0.])
        lr_err = array.array('d', [0.])
        tree.SetBranchAddress(branch_lr,     lr)
        tree.SetBranchAddress(branch_lr_err, lr_err)
        nof_entries = tree.GetEntries()
        for i in range(nof_entries):
          tree.GetEntry(i)
          sig_bkg_dict['lr'].append(lr[0])
          sig_bkg_dict['lr_up'].append(lr[0] + lr_err[0])
          sig_bkg_dict['lr_down'].append(lr[0] - lr_err[0])
      else:
        branch_sig_err = '%s_err' % branch_sig
        branch_bkg_err = '%s_err' % branch_bkg
        missing_branches = {branch_sig, branch_sig_err, branch_bkg, branch_bkg_err} - branch_names
        if missing_branches:
          raise ValueError("No such branches in {tree_name}: {missing_branches}".format(
            tree_name=tree_name,
            missing_branches=' '.join(missing_branches)
          ))

        # read the branches, assuming type double
        sig     = array.array('d', [0.])
        sig_err = array.array('d', [0.])
        bkg     = array.array('d', [0.])
        bkg_err = array.array('d', [0.])
        tree.SetBranchAddress(branch_sig,     sig)
        tree.SetBranchAddress(branch_sig_err, sig_err)
        tree.SetBranchAddress(branch_bkg,     bkg)
        tree.SetBranchAddress(branch_bkg_err, bkg_err)
        nof_entries = tree.GetEntries()
        for i in range(nof_entries):
          tree.GetEntry(i)
          lr, lr_err = 0., 0.
          if (sig[0] + bkg[0]) > 0.:
            lr     = sig[0] / (sig[0] + bkg[0])
            lr_err = np.sqrt((sig[0] * sig_err[0])**2 + (bkg[0] * bkg_err[0])**2) / (sig[0] + bkg[0])**2
          sig_bkg_dict['lr'].append(lr)
          sig_bkg_dict['lr_up'].append(lr + lr_err)
          sig_bkg_dict['lr_down'].append(lr - lr_err)

      # now that even though the TTree likely contains the up/down values of the likelihood ratio
      # these are kind-of systematics which do not reflect the accuracy of the inegration very well
      # therefore we compute the up-down values ourselves
      sig_bkg_dict['lr']      = np.asarray(sig_bkg_dict['lr'])
      sig_bkg_dict['lr_up']   = np.asarray(sig_bkg_dict['lr_up'])
      sig_bkg_dict['lr_down'] = np.asarray(sig_bkg_dict['lr_down'])

  logging.debug("Processing the input data")
  # now we create three set of ROC curves for each signal-background pair
  # 1) using lr for both signal and background
  # 2) using lr + err for signal and lr - err for background (to get lower bound)
  # 3) using lr - err for signal and lr + err for background (to get upper bound)
  # the output file names contain signal and background sample names, and the branch names used to create the LR

  nof_wps = 100
  wps = np.linspace(0., 1., nof_wps + 1)

  for label in input:
    for sig_bkg in input[label]:
      sig_bkg['default'] = get_roc(wps, sig_bkg['lr'])
      sig_bkg['up']      = get_roc(wps, sig_bkg['lr_up'])
      sig_bkg['down']    = get_roc(wps, sig_bkg['lr_down'])

  # merge the signal efficiencies of each version of likelihood ratio
  # compute missing background effiencies by linearly interpolating b/w nearest background efficiencies
  sig = input['sig'][0]
  sig_effs = np.sort(np.unique(np.concatenate([sig['default'], sig['down'], sig['up']])))
  for bkg in input['bkg']:
    bkg['default'] = np.interp(sig_effs, sig['default'][::-1], bkg['default'][::-1])[::-1]
    bkg['down']    = np.interp(sig_effs, sig['up'][::-1],      bkg['down'][::-1])[::-1]
    bkg['up']      = np.interp(sig_effs, sig['down'][::-1],    bkg['up'][::-1])[::-1]
    bkg['wp_default'] = np.interp(sig_effs, sig['default'][::-1], wps[::-1])[::-1]
    bkg['wp_up']      = np.interp(sig_effs, sig['down'][::-1],    wps[::-1])[::-1]
    bkg['wp_down']    = np.interp(sig_effs, sig['up'][::-1],      wps[::-1])[::-1]
  sig_effs = sig_effs[::-1]

  # save the results
  for bkg in input['bkg']:
    output_basename = 'roc_{signal}_{background}_{branch_name}.csv'.format(
      signal      = sig['unique_id'],
      background  = bkg['unique_id'],
      branch_name = branch_lr if branch_lr else '-'.join([branch_sig, branch_bkg]),
    )
    output_filename = os.path.join(output_dir, output_basename)
    np.savetxt(
      output_filename,
      np.transpose(np.asarray([
        sig_effs, bkg['default'], bkg['wp_default'],
                  bkg['up'],      bkg['wp_up'],
                  bkg['down'],    bkg['wp_down'],
      ])),
      fmt       = '%.5f',
      delimiter = ',',
    )
    logging.debug("Created file {filename}".format(filename = output_filename))
