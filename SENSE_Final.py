#!/usr/bin/python
# SENSE_Final.py -s /home/mcarias/Marc_MRI/Cylindrical/LiveJuly5/3D_LFOVlive_2013_07_05.fid -b 3D_SFOVlive_2013_07_25.fid -o foo.mnc

from sys import path as syspath
syspath.append('/home/bjnieman/source/vnmr')
from varian_recon import *
from varian_read_file import *
from pylab import*
from scipy.interpolate import splrep,splev
import optparse

class dummy:
    def __init__(self,r,i,p,m,f2d,f3d,down_RO,max_set,range_max,range_min,large_data,shiftro,shiftpe,shiftpe2):
        self.vType=None
        self.real=r
        self.imag=i
        self.phase=p
        self.mag=m
        self.fft2d=f2d
        self.fft3d=f3d
        self.downsample_RO=down_RO
	self.max_range=max_set
        self.image_range_max=range_max
        self.image_range_min=range_min
        self.large_data_recon=large_data
        self.fov_shift_ro=shiftro
        self.fov_shift_pe1=shiftpe
        self.fov_shift_pe2=shiftpe2

def phase_corr_w0reps(kdata,petable_name,petable_arrays=('t1','t2')):
    enc_array_1 = parse_petable_file(petable_name,petable_arrays[0])
    enc_array_1 = enc_array_1
    enc_array_2 = parse_petable_file(petable_name,petable_arrays[1])
    enc_array_2 = enc_array_2
    i1,i2=nonzero(reshape((enc_array_1==0)*(enc_array_2==0),(len(enc_array_1)/kdata.shape[1],kdata.shape[1])))
    reps = prod(kdata.shape[0:2])/len(enc_array_1)
    if (reps>1):
        for q in range(1,reps):
            i1=append(i1,i1[0:len(i1)/q]+q*len(enc_array_1)/kdata.shape[1])
            i2=append(i2,i2[0:len(i2)/q])
    zerfids = kdata[i1,i2,:]
    maxind = argmax(abs(zerfids[0,:]))
    #spline fit and phase determination
    kwidth=(max(i1)*kdata.shape[1]+max(i2))/5
    kpts=arange(0,kwidth,4.5*kwidth)+kwidth/2
    splpar=splrep(i1*kdata.shape[1]+i2,zerfids[:,maxind].real,t=kpts)
    splfitreal=splev(arange(kdata.shape[0]*kdata.shape[1]),splpar)
    splpar=splrep(i1*kdata.shape[1]+i2,zerfids[:,maxind].imag,t=kpts)
    splfitimag=splev(arange(kdata.shape[0]*kdata.shape[1]),splpar)
    phval=arctan2(splfitimag,splfitreal)
    #phase correction
    kdata=kdata*exp(-1.j*reshape(phval,(kdata.shape[0],kdata.shape[1])))[:,:,newaxis]
    return kdata

def main():
    usage = "%prog [options] --input-dir INPUTDIR"
    p = optparse.OptionParser(usage=usage)
    p.add_option('--input-dir-lr', '-s', dest="indirlr", help="input directory")
    p.add_option('--input-dir-hr', '-b', dest="indirhr", help="input directory")
    p.add_option('--outfile'  , '-o', dest="outfile", help="output")
    ops, arguments = p.parse_args()

    inputdirectory=ops.indirlr

    petable='/home/mcarias/Marc_MRI/table_test/75_75repk'
    vnmrfidfilelist,data_shape,header_info,param_dict,procpar_text_lines = open_vnmrfid_file(inputdirectory)
    close_vnmrfid_file(vnmrfidfilelist)
    nmice=7 #param_dict['nmice']
    #image_data=zeros([nmice,param_dict['nv'],param_dict['nv2'],(param_dict['np']/2)],complex)
    image_data=zeros([nmice,param_dict['nv2'],param_dict['nv'],(param_dict['np']/2)],complex)
    
    for j in range(nmice):
        temp,param_dict,procpar_text_lines = gen_kspace_simple(inputdirectory,j)
        temp = petable_orderedpair_reordering(temp,petable_arrays=('t1','t2'),petable_name=petable,matrix=(param_dict['nv'],param_dict['nv2']))
        image_data[j,:,:,:] = temp.copy()
    
    # shift in z
    ro_start=-10.0
    ro_end=10.0
    ro_step=0.5
    ro_pos=arange(ro_start,ro_end,ro_step)
    ropixshift=[0]
    ro_pixshift = zeros((7,),float) 
    for j in [0,1]:
        img0 = recon_3d(image_data[j,:,:,:])
        for k in [[2,5],[3,4,6]][j]:
            Cval=[]
            thresh = nan
            for q in ro_pos:
                temp=fov_shift(image_data[k,:,:,:],float(q),-1,False)
                temp=recon_3d(temp)
                if isnan(thresh):
                    thresh = sort(ravel(abs(temp)))[-int(0.01*product(temp.shape))]
                i1,i2,i3 = nonzero(abs(temp)>thresh)
                mask = nonzero(abs(temp)>thresh)
                Nmask = sum(mask)
                Cval.append(abs(sum((img0[i1,i2,i3]-sum(img0[i1,i2,i3])/Nmask)*conj(temp[i1,i2,i3]-sum(temp[i1,i2,i3])/Nmask))))
            Cval=array(Cval)
            ind = argmax(Cval)
            p = polyfit(ro_pos[ind-3:ind+3],Cval[ind-3:ind+3],2)
            ro_pixshift[k] = -0.5*p[1]/p[0]



#######################################
#try deriving only from SFOV image
#######################################

    inputdirectory=ops.indirhr

    petable='/micehome/mcarias/Marc_MRI/table_test/222_222rep3'
    vnmrfidfilelist,data_shape,header_info,param_dict,procpar_text_lines = open_vnmrfid_file(inputdirectory)
    close_vnmrfid_file(vnmrfidfilelist)
    nmice=7 #param_dict['nmice']
    image_data=zeros([nmice,param_dict['nv2'],param_dict['nv'],(param_dict['np']/2)],complex)
    for j in range(nmice):
        print "Mouse %d / %d ..."%(j,nmice)
        temp,param_dict,procpar_text_lines = gen_kspace_simple(inputdirectory,j)
        temp = phase_corr_w0reps(temp,petable)
        temp = petable_orderedpair_reordering(temp,petable_arrays=('t1','t2'),petable_name=petable,matrix=(param_dict['nv'],param_dict['nv2']))
        #temp = fov_shift(temp,ro_pixshift[j],-1,0)
        image_data[j,:,:,:] = fermi_ellipse_filter(temp)

#    C=ones((nmice,nmice),complex)
#    pe1ramp=(arange(data_shape[-2])-data_shape[-2]/2).astype(float)/(data_shape[-2]/2)
#    pe2ramp=(arange(data_shape[-3])-data_shape[-3]/2).astype(float)/(data_shape[-3]/2)
#    I1pe1 = zeros(image_data.shape[1:],image_data.dtype)
#    I1pe2 = zeros(image_data.shape[1:],image_data.dtype)
#    I2pe1 = zeros(image_data.shape[1:],image_data.dtype)
#    I2pe2 = zeros(image_data.shape[1:],image_data.dtype)
#    for j in range(nmice):
#        for k in range(j+1,nmice):
#            I1pe1 = recon_3d(image_data[j,:,:,:]*pe1ramp[newaxis,:,newaxis])
#            I2pe1 = recon_3d(image_data[k,:,:,:]*pe1ramp[newaxis,:,newaxis])
#            I1pe2 = recon_3d(image_data[j,:,:,:]*pe2ramp[:,newaxis,newaxis])
#            I2pe2 = recon_3d(image_data[k,:,:,:]*pe2ramp[:,newaxis,newaxis])
#            abs1 = sqrt(abs(I1pe1)**2+abs(I1pe2)**2)
#            abs2 = sqrt(abs(I2pe1)**2+abs(I2pe2)**2)
#            mask1 = abs1 > sort(ravel(abs1))[-int(0.008*prod(image_data.shape[1:]))]
#            mask2 = abs2 > sort(ravel(abs2))[-int(0.008*prod(image_data.shape[1:]))]
#            C[j,k] = 0.5*sum(I1pe1*conj(I2pe1)*mask2)/sum(I2pe1*conj(I2pe1)*mask2) + 0.5*sum(I1pe2*conj(I2pe2)*mask2)/sum(I2pe2*conj(I2pe2)*mask2)
#            C[k,j] = 0.5*sum(I2pe1*conj(I1pe1)*mask1)/sum(I1pe1*conj(I1pe1)*mask1) + 0.5*sum(I2pe2*conj(I1pe2)*mask1)/sum(I1pe2*conj(I1pe2)*mask1)
#    
#    from scipy.linalg import inv
#    C_inv = inv(C)
#
#
#
#
#    imgrecon = zeros(image_data.shape,image_data.dtype)
#    for k in range(nmice):
#        cimg =(image_data[k,:,:,:])
#        for j in range(nmice):
#            imgrecon[j,:,:,:] += C_inv[j,k]*cimg
#    
    
#    shift=param_dict['mm_ppe']*param_dict['nv']/param_dict['lpe']
#    shift2=param_dict['mm_ppe2']*param_dict['nv2']/param_dict['lpe2']
#    shift3=[0.0,0.0,100.0,0.0,0.0,0.0,0.0]
#    for j in range(7):
#        imgrecon[j,:,:,:] = fov_shift(imgrecon[j,:,:,:],shift2[j],-3,0)
#    
#    for j in range(7): 
#        imgrecon[j,:,:,:] = fov_shift(imgrecon[j,:,:,:],shift[j],-2,0)
#    
#    for j in range(7): 
#        imgrecon[j,:,:,:] = fov_shift(imgrecon[j,:,:,:],shift3[j],-1,0)
    
    
#    for j in range(7): 
#        imgrecon[j,:,:,:] = recon_3d(imgrecon[j,:,:,:]).astype(N.complex)


    for j in range(7): 
        image_data[j,:,:,:] = recon_3d(image_data[j,:,:,:]).astype(N.complex)



    options=dummy(0,0,0,1,0,1,0,0,0,0,0,0.0,0.0,0.0)
    write_to_mnc_file(ops.outfile,image_data,param_dict,procpar_text_lines,options)

#array([  0.        ,   0.        ,  -3.15715602,  15.27943306,
#         5.48314156,  -6.92231147,  -6.99088097]) #july25
#
#array([ 0.        ,  0.        ,  0.91471226, -0.83587943,  9.50157529,
#       -3.49625152, -3.25224687]) #july5
#
#
#array([ 0.        ,  0.        ,  6.08231155, -1.12969368,  7.33308884,
#        1.4385015 , -3.44495199])
#
#
#
#
#
#
#

if __name__ == '__main__':
    main()
    sys.exit(0)

