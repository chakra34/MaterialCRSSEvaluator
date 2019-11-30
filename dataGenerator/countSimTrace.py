#!/usr/bin/env python

# ----- ABBREVIATION ----- #
# ss: Slip System
# sf: Schmid Factor
# gamma: accumulated shear
# r: ratio
# idx/idc: index/indices
# vctr: vector
# gid: grain ID

import os
import damask
import numpy as np
from optparse import OptionParser

scriptName = os.path.splitext(os.path.basename(__file__))[0]
scriptID = ' '.join([scriptName, damask.version])


# ----- Helper functions ----- #
def slipTraceIsVisiblility(gammas_ndarray, threshold, weights=None):
    """Return visibility vector for given grain."""
    # NOTE: if the accumulated shear for given slip system accounts for
    #       over 20% of the total shear, it is considered to have a visible
    #       slip trances in this grain.
    sum_gamma = sum(map(abs, gammas_ndarray))
    r_gamma = abs(gammas_ndarray)/sum_gamma

    # weight the threshold by given weights if possible
    if weights is not None:
        threshold = threshold * np.array(weights)

    return r_gamma > threshold


# ----- MAIN ----- #
desp_msg = "Generate observation&population list from CPFFT simulation."
parser = OptionParser(option_class=damask.extendableOption,
                      usage='%prog options [file[s]]',
                      description=desp_msg,
                      version=scriptID)
parser.add_option('-t', '--threshold',
                  dest='threshold',
                  type=float,
                  help='threshold for slip trace count')
parser.add_option('-w', '--weighted',
                  dest='weighted',
                  action='store_true',
                  help='weight the threshold by slip family size')

# this is a modified script, always calculate all possible thresholds
parser.set_defaults(threshold=0.2,
                    weighted=False)

(options, filenames) = parser.parse_args()

# ----- going through list of files ----- #
for fname in filenames:
    # ----- read in data using DAMASK ASCII table class ----- #
    asciiTable = damask.ASCIItable(name=fname, buffered=True)
    asciiTable.head_read()
    asciiTable.data_readArray()  # really slow for large dataset

    # getting texture data
    texture_idx = asciiTable.label_index('texture')
    texture_vctr = np.copy(asciiTable.data[:, texture_idx])

    # getting accumulated shear data
    gamma_idx = asciiTable.label_index('accumulatedshear_slip')
    gamma_dim = asciiTable.label_dimension('accumulatedshear_slip')
    gamma_idc = np.arange(gamma_idx, gamma_idx+gamma_dim)
    gammas = np.copy(asciiTable.data[:, gamma_idc])

    # getting orientations data
    qs_idx = asciiTable.label_index('orientation')
    qs_dim = asciiTable.label_dimension('orientation')
    qs_idc = np.arange(qs_idx, qs_idx+qs_dim)
    qs = np.copy(asciiTable.data[:, qs_idc])

    # getting the Schmid factor
    # NOTE: Unfortunately Schmid factor is not recognized as a group
    #       quantity.
    # NOTE: CPFFT simulation does not have second order prism slip and
    #       the last 6 pyramidal <c+a> slip.
    ss_idc_selected = np.array(range(6) + range(9, 27))
    for label in asciiTable.labels():
        if "S[" in label:
            if asciiTable.label_dimension(label) > 1:
                sf_idx = asciiTable.label_index(label)
                sf_dim = asciiTable.label_dimension(label)
                sf_idc = np.arange(sf_idx, sf_idx+sf_dim)[ss_idc_selected]
                break
    sfs = np.copy(asciiTable.data[:, sf_idc])

    # close table
    asciiTable.close(dismiss=True)

    # ----- count slip traces in each grain ----- #
    gids = set(texture_vctr)  # get unique grain id
    gammas_gid = []
    sfs_gid = []
    qs_gid = []
    # calculate grain average accumulated shear and Schmid factors
    for gid in gids:
        gid_idc = np.where(texture_vctr == gid)[0]
        gammas_gid.append(np.mean(gammas[gid_idc, :], axis=0))
        sfs_gid.append(np.mean(sfs[gid_idc, :], axis=0))
        # ----- get grain average orientation ----- #
        qs_tmp = qs[gid_idc, :]
        qs_tmp = [damask.Orientation(quaternion=damask.Quaternion(quatArray=q))
                  for q in qs_tmp]
        qs_gid.append(damask.Orientation.average(qs_tmp).asQuaternion())

    gammaSum_gid = [sum(map(abs, gamma)) for gamma in gammas_gid]
    with open(fname.replace(".txt", ".grainSum.txt"), 'w') as f:
        outstr = '1  header\n gammaSum\n'
        outstr += "\n".join(map(str, gammaSum_gid))
        f.write(outstr)

    # slip family ID
    sf_ids = [1]*3 + [2]*3 + [3]*6 + [4]*12

    # write out the population file
    # --> all slip systems in all grains
    with open(fname.replace(".txt", ".population.txt"), 'w') as f:
        outstr = '1  header\nFamily\tSF\t'
        outstr += '\t'.join(["{}_q".format(i+1) for i in xrange(4)]) + '\n'

        for i, gid in enumerate(gids):
            q_str = "\t".join(map(str, qs_gid[i]))
            outstr += "\n".join(['{}\t{}\t'.format(fid, sf) + q_str
                                for fid, sf in zip(sf_ids,
                                                   sfs_gid[i])])
            outstr += '\n'
        f.write(outstr)

    # count slip trace for each grain
    # --> only slip systems that satisfy slip trace condition
    options.threshold = [0.01, 0.1, 0.2, 0.3, 0.4]
    # 1 for basal/prism, 0.5 for pyra, 0.25 for pyrca
    if options.weighted:
        weights = np.array([1]*6 + [0.5]*6 + [0.25]*12)
    else:
        weights = None
    for threshold in options.threshold:
        outfname = fname.replace(".txt", ".observation")
        outfname += "_" + str(threshold) + ".txt"
        with open(outfname, 'w') as f:
            outstr = '1  header\nFamily\tSF\t'
            outstr += '\t'.join(["{}_q".format(i+1) for i in xrange(4)]) + '\n'
            for i, gid in enumerate(gids):
                visible = slipTraceIsVisiblility(gammas_gid[i],
                                                 threshold,
                                                 weights=weights)
                q_str = "\t".join(map(str, qs_gid[i]))
                if any(visible):
                    outstr += "\n".join(["{}\t{}\t".format(sf_ids[j],
                                                           sfs_gid[i][j])+q_str
                                         for j, viz in enumerate(visible)
                                         if viz])
                    outstr += '\n'
            f.write(outstr)
