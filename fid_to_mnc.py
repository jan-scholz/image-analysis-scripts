# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 13:17:26 2015
Script to reconstruct 2D Varian data (which won't reconstruct in vrecon)
@author: sportnoy
"""
import numpy as np
import re
import os
import sys
from pyminc.volumes.factory import volumeFromDescription
from optparse import OptionParser



def procpar_read(parfile):
    """Function to read in procpar file and create a parameter dictionary"""
    f = open(parfile,'r')   
    procpar_text_lines = f.readlines()
    f.close()
    param_dict = procpar_text_lines_to_dict(procpar_text_lines)
    return param_dict


def readraw(foldername):
    """Function to parse varian fid data"""
    datafile=foldername+'/fid'
    parfile=foldername+'/procpar'
    if options.outfilename:
        imgfile=options.outfilename
    else:    
        imgfile=foldername+'/img.mnc'
    phaseimgfile=foldername+'/phase.mnc'
    param_dict=procpar_read(parfile)
    
    seqcon=param_dict['seqcon']
    nslices=param_dict['ns']
    
    f=open(datafile,'rb') 
    nblocks = np.fromfile(f,np.dtype('>I4'),1) #number of blocks in file, read BIG-ENDIAN
    ntraces = np.fromfile(f,np.dtype('>I4'),1)
    nelems  = np.fromfile(f,np.dtype('>I4'),1)
    ebytes  = np.fromfile(f,np.dtype('>I4'),1) #number of bytes per element
    tbytes  = np.fromfile(f,np.dtype('>I4'),1) #number of bytes per trace
    bbytes  = np.fromfile(f,np.dtype('>I4'),1) #number of bytes per block

    f.seek(8,1)
       
    kspc = np.zeros([nblocks,ntraces,nelems/2],dtype = 'complex')
    img  = np.zeros([nblocks,ntraces,nelems/2],dtype = 'complex')
    count=0
    while count<nblocks:
        f.seek(28,1)
        raw_data = np.fromfile(f,np.dtype('>f4'),nelems*ntraces,'')  
        complex_data = raw_data[::2]+1j*raw_data[1::2]
        complex_data.shape = (ntraces,nelems/2)
        kspc[count,:,:] = complex_data
        img[count,:,:] = np.fft.fftshift(np.fft.ifft2(complex_data))
        count = count+1
    f.close()   
    
    #Parse data with self-gating
    #if(get_dict_value(param_dict,'sgflag','n')=='y'):   
    if(param_dict.get('sgflag','n')=='y'):    
        print 'self gating'
        sg_fids,kspc_new = get_sg_signal(kspc,param_dict)
        img = np.zeros(kspc_new.shape,dtype = 'complex')
        count = 0
        while count<nblocks:
            img[count,:,:] = np.fft.fftshift(np.fft.ifft2(kspc_new[count,:,:]))
            count = count+1
    
    #Parse multi-slice data
    if(seqcon=='nccnn' and nslices>1):
    #need to reparse k-space data and ifft, as all slices will be stored in a single block
       if nblocks>1:      
           block_count=0
           kspc_new=np.zeros([nblocks,nslices,ntraces/nslices,nelems/2],dtype='complex')
           img_new=np.zeros([nblocks,nslices,ntraces/nslices,nelems/2],dtype='complex')
           while block_count<nblocks:
                slice_count=0
                while slice_count<nslices:
                    kspc_new[block_count,slice_count,:,:]=np.squeeze(kspc)[block_count,slice_count::nslices,:]
                    img_new[block_count,slice_count,:,:]=np.fft.fftshift(np.fft.ifft2(kspc_new[block_count,slice_count,:,:]))
                    slice_count=slice_count+1
         
                block_count=block_count+1
                
       else:
           kspc_new=np.zeros([nslices,ntraces/nslices,nelems/2],dtype='complex')
           img_new=np.zeros([nslices,ntraces/nslices,nelems/2],dtype='complex')
           slice_count=0
           while slice_count<nslices:
                kspc_new[slice_count,:,:]=np.squeeze(kspc)[slice_count::nslices,:]
                img_new[slice_count,:,:]=np.fft.fftshift(np.fft.ifft2(kspc_new[slice_count,:,:]))
                slice_count=slice_count+1
       kspc,img=kspc_new,img_new         
    
    #Parse 2D fourier flow encoded data 
    #if(get_dict_value(param_dict,'fourier_flow','n')=='y'): 
    if(param_dict.get('fourier_flow','n')=='y'):    
       #nve=get_dict_value(param_dict,'num_velocity_encodes',0)
       nve=param_dict.get('num_velocity_encodes',0)
       if(nve==0):
           #nve=get_dict_value(param_dict,'numarray',1)
           nve=param_dict.get('numarray',1)
       print nve
       kspc_new=np.zeros([nve,ntraces/nve,nelems/2],dtype='complex')
       img_new=np.zeros([nve,ntraces/nve,nelems/2],dtype='complex')
       velocity_count=0
       #if(get_dict_value(param_dict,'fourier_loopmode','inner')=='inner'):
       if(param_dict.get('fourier_loopmode','inner')=='inner'):   
           while velocity_count<nve:
               kspc_new[velocity_count,:,:]=np.squeeze(kspc)[velocity_count::nve,:]
               img_new[velocity_count,:,:]=np.fft.fftshift(np.fft.ifft2(kspc_new[velocity_count,:,:]))
               velocity_count=velocity_count+1 
       else:
          while velocity_count<nve:
               kspc_new[velocity_count,:,:]=np.squeeze(kspc)[velocity_count*ntraces/nve:(velocity_count+1)*ntraces/nve,:]
               img_new[velocity_count,:,:]=np.fft.fftshift(np.fft.ifft2(kspc_new[velocity_count,:,:]))
               velocity_count=velocity_count+1  
               
       kspc,img=kspc_new,img_new
    
    #Parse 2D velocity phase contrast data   
    #if(get_dict_value(param_dict,'fc','n')=='y' and not get_dict_value(param_dict,'vencz',0)==0):
    if(param_dict.get('fc','n')=='y') and not (param_dict.get('vencz',0)==0):
        print 'velocity phase contrast'        
        kspc_new=np.zeros([2,ntraces/2,nelems/2],dtype='complex')
        img_new=np.zeros([2,ntraces/2,nelems/2],dtype='complex')
        count=0
        while count<2:
            kspc_new[count,:,:]=np.squeeze(kspc)[count::2]
            img_new[count,:,:]=np.fft.fftshift(np.fft.ifft2(kspc_new[count,:,:]))
            count=count+1
            
        kspc,img=kspc_new,img_new
    
    if(options.magnitude):           
        write_to_mnc_file_coords(imgfile,np.abs(img),param_dict,'mag')

    if(options.phase):
        write_to_mnc_file_coords(phaseimgfile,np.abs(img),param_dict,'phase')             

    return param_dict,np.squeeze(kspc),np.squeeze(img)
    

def write_to_mnc_file_coords(filename,data,param_dict,imgtype='mag'):
    """Function to write to a minc file with coordinates"""
    out_im=data    
    vType='ushort'
    if(imgtype=='phase'):
        vType='short'
    dim_sizes=np.shape(out_im)
    nD=len(dim_sizes)
    fov_array = [10.0*float(param_dict.get(x,1.0)) for x in ['lpe2','lpe','lro']]
    fov_ro = fov_array[-1]
    fov_pe = fov_array[-2]
    fov_pe2 = fov_array[-3]    
    ns=param_dict.get('pss',0.0).size
    

    if param_dict.get('nD',1.0)==2.0:    
	#loop order: slice,pe1,ro

        #get offset information:
        mppe=10.0*param_dict.get('mppe',0.0)[0]
        mpss=10.0*param_dict.get('mpss',0.0)[0]
        mpro=10.0*param_dict.get('mpro',0.0)[0]
        
        if (len(np.shape(np.squeeze(out_im)))==2 
                                or ns==1):
            pss=10.0*param_dict.get('pss',0.0)
            start_z=pss+mpss
            step_z=param_dict.get('thk',0.0)
             
        else:
            pss=10.0*param_dict.get('pss',np.zeros(dim_sizes[-3]))
            start_z=pss[0]+mpss
            step_z=pss[1]-pss[0]
          
        step_y=fov_pe/dim_sizes[-2]
      
        step_x=fov_ro/dim_sizes[-1]            
            
        oepf=param_dict.get('oepf','n')
        
        if(oepf=='y'):
            start_y=-mppe-fov_pe/2.0 
        else:
            start_y=mppe-fov_pe/2.0

        start_x=mpro-fov_ro/2.0
        
    else:
        start_z=-fov_pe2/2.0
        start_y=-fov_pe/2.0
        start_x=-fov_ro/2.0
        
        step_z=fov_pe2/dim_sizes[-3]
        step_y=fov_pe/dim_sizes[-2]
        step_x=fov_ro/dim_sizes[-1]

    start_tup=(start_z,start_y,start_x,0,0)
    step_tup=(step_z,step_y,step_x,1,1)
    
    start_tup=(0,0,0,0,0)    
    step_tup=(1,1,1,1,1)
    
    orient=param_dict.get('orient','trans')
    if (orient=='trans' or  orient=='oblique'):
        dim_names=tuple(['zspace','xspace','yspace'])
        start_tup=(start_z,-1.0*start_y,start_x,0,0)
        step_tup=(step_z,-1.0*step_y,step_x,0,0)
        print 'trans'
        
    if (orient=='sag'):
        dim_names=tuple(['xspace','yspace','zspace'])
        start_tup=(-1.0*start_z,start_y,start_x,0,0)
        step_tup=(-1.0*step_z,step_y,step_x,0,0) 
        print 'sag'
    
    if (orient=='cor'):
        dim_names=tuple(['yspace','xspace','zspace'])
        start_tup=(-1.0*start_z,-1.0*start_y,start_x,0,0)
        step_tup=(-1.0*step_z,-1.0*step_y,step_x,0,0) 
        print 'cor'
   
   
    if (nD<4):
        nout_files=1
    else:
       nout_files=dim_sizes[0]
       nD = 3

    for j in range(nout_files):
        if (nout_files==1):
            outfile_name = filename
            imgshape = dim_sizes
        else:
            outfile_name = filename[:-4] + '_' + str(j) + '.mnc'
            imgshape = dim_sizes[1:]
        print "Writing to output file %s..." % outfile_name

        if os.path.exists(outfile_name):
            os.remove(outfile_name) 
        
        mnc_vol = volumeFromDescription(outfile_name,dim_names,imgshape,start_tup[0:nD],
                                            step_tup[0:nD],volumeType=vType,dtype=out_im.dtype)
        if (nout_files==1):
            mnc_vol.setdata(out_im)
        else:
            mnc_vol.setdata(out_im[j])
        mnc_vol.writeFile()
        mnc_vol.closeVolume() 
    return nout_files 


def procpar_text_lines_to_dict(text_lines):
    """Function to create a parameter dictionary from the text in the procpar file.
       Obtained from Brian Nieman"""
    param_dict={}
    ntext = 0
    curr_line = 0
    nlines = len(text_lines)
    while (curr_line<nlines):
        line = text_lines[curr_line]
        m = re.search('^[a-z]',line)
        if not m:
            curr_line += 1
            continue
        else:
            words = line.split()
            varname = words[0]
            subtype = int(words[1])
            basictype = int(words[2])
            curr_line += 1
            line = text_lines[curr_line]
            words = line.split()
            nvals = int(words[0])
            if (basictype==1) and (subtype==7) and (not varname=='filter') and (not varname=='dres'):
                try:
                    varval = np.array([int(x) for x in words[1::]],int)
                except ValueError:
                    varval = np.array([float(x) for x in words[1::]],float)
            elif (basictype==1):
                varval = np.array([float(x) for x in words[1::]],float)
            else:
                varval = [ re.search('".*"',line).group()[1:-1] ]
                while (len(varval)<nvals):
                    curr_line+=1
                    line = text_lines[curr_line]
                    varval.append( re.search('".*"',line).group()[1:-1] )
            if (nvals==1): varval=varval[0]
            param_dict[varname] = varval
            curr_line+=1
    return param_dict          

if __name__ == '__main__':
#    foldername = sys.argv[-1]
#    readraw(foldername)    
    parser=OptionParser()
    parser.add_option("--mag", default=1, action="store_true", dest="magnitude", 
                      help="reconstruct magnitude image")
    parser.add_option("--phase", default=0, action="store_true", dest="phase",
                      help="reconstruct phase image") 
    parser.add_option("--outputfile", dest="outfilename", 
                      help="reconstructed minc filename")                  
    options,args=parser.parse_args()
    foldername=args[-1]
    readraw(foldername)
