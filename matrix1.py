#!/usr/bin/env python

import sys
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt

import resource
rusage_denom = 1024.

# hexagonal subdivision of sphere
# calc area in white surface space
# calc mean area
# calac deviation from mean area
# adjust connectivity by area of each piece?
# adjust by sulcal depth


# connected_components(G)[0]  # get biggest connected component


# 29078 30069
# sed -n '29078,35083p' fdt_matrix1.dot > foo.dot
# mris_info S08/surf/lh.white  # nvertices   122506
# unique in dot 122506
# total 49506200
#
# mri_surf2surf --hemi lh --srcsubject 'subjectID' --sval-xyz pial 
# --trgsubject fsaverage6 --trgicoorder 6 --trgsurfval 'outputfilename'
# http://webcache.googleusercontent.com/search?q=cache:bNe3c_iU0_MJ:https://www.jiscmail.ac.uk/cgi-bin/webadmin%3FA3%3Dind1302%26L%3DFSL%26E%3Dquoted-printable%26P%3D6109327%26B%3D--Apple-Mail-299-911161879%26T%3Dtext%252Fhtml%3B%2520charset%3Dus-ascii%26attachment%3Dq%26XSS%3D3+&cd=10&hl=en&ct=clnk
# https://www.jiscmail.ac.uk/cgi-bin/webadmin?A3=ind1207&L=FSL&E=quoted-printable&P=1274004&B=--_f6d846f2-39ed-4c5f-bfb6-abeffcc68f89_&T=text%2Fhtml;%20charset=iso-8859-1&pending=
#
# creates surfaces for subcortical structures
# aseg2srf -s subject
# http://brainder.org/2012/05/08/importing-freesurfer-subcortical-structures-into-blender/


if __name__ == "__main__":
    m = sys.argv[1]
    X = np.loadtxt(m);
    print "X shape", X.shape
    print type(X)
    #print X


    # convert spartse (where) matrix to dense matrix
    sB = sparse.coo_matrix((X[:,2],(X[:,0],X[:,1])))
    print 'sparse, mem:', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.

    B = np.array(sB.todense())
    print 'dense mem:', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.

    B = B[0:min(B.shape),0:min(B.shape)]
    print "B shape", B.shape
    print type(B)
    #print B

    fig = plt.figure()
    ax = plt.subplot(111)
    #ax.axis('equal')

    ax.imshow(B,interpolation='nearest')
    #plt.show()
    #plt.show(block=True)
    fig.savefig('plot.png',dpi=100)


