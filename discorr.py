#!/usr/bin/python
# formally known as distortion_correction_november2007.pl
# modified by jan.scholz@phenogenomics.ca

# Takes a brain MR scan and will apply the november 2007 distortion
# correction transformation to it:
#
# -MR to CT transform
#
# Will determine the coil the mouse was scanned in 
# using the vnmr:coil entry in the MINC header

# TODO:
#   - return queue number
#   - job array

version = "0.1";

import sys
import os
import optparse
import subprocess
#import argparse

#/projects/mice/share/arch/linux-x86_64-eglibc2_11_1/bin/sge_batch -q bigmem.q,all.q -l vf=4G -J dc-4G-rot36_b0 mincresample -sinc -2 -transform /projects/mice/matthijs/distortion-correction-november-2007/coil3-phantom-16oct2007/coil3-mri-to-CAD-lsq6-nlin-extra-large-super-inverted-lsq6.xfm -like ../raw/rot03_b0.mnc ../raw/rot36_b0.mnc discorr/rot36_b0.november_2007_distortion_corrected.mnc





basedir = "/projects/mice/matthijs/distortion-correction-november-2007/"

coil1_xfm = basedir + "/coil1-phantom-17oct2007/coil1-mri-to-CAD-lsq6-nlin-extra-large-super-inverted-lsq6.xfm"
coil2_xfm = basedir + "/coil2-phantom-18oct2007/coil2-mri-to-CAD-lsq6-nlin-extra-large-super-inverted-lsq6.xfm"
coil3_xfm = basedir + "/coil3-phantom-16oct2007/coil3-mri-to-CAD-lsq6-nlin-extra-large-super-inverted-lsq6.xfm"

coil1_likefile = basedir + "/coil1-phantom-17oct2007/img_17oct07.0-reshaped-coil1-likefile.mnc"
coil2_likefile = basedir + "/coil2-phantom-18oct2007/img_18oct07.1-reshaped-coil2-likefile.mnc"
coil3_likefile = basedir + "/coil3-phantom-16oct2007/img_16oct07.2-reshaped-coil3-likefile.mnc"


def validcoil(x):
	try:
		coil = int(x)
	except ValueError:
		raise
	
	if not coil in range(1,17):
		raise ValueError

	return coil


def splitarg(a):
	filename = a.split(":")[0]
	coil = float("nan")

	try:
		coil = validcoil(a.split(":")[1])
	except IndexError:
		pass
	except ValueError:
		print >> sys.stderr, 'ERROR: not a valid coil number:', a.split(":")[1]
		sys.exit(1)

	return (filename, coil)
	

def getcoil(filename):
	try:
		p = subprocess.Popen(['mincinfo','-attvalue','vnmr:coil',filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	except OSError:
		pass

	if p.wait():
		if p.communicate()[1].startswith('Error reading file.'):
			coil = float("nan")
		else:
			print >> sys.stderr, 'ERROR: could not read file', filename
			sys.exit(1)
	else:
		try:
			coil = int(p.communicate()[0].strip())
		except ValueError:
			print >> sys.stderr, 'ERROR: read invalid coil number from file:', filename
			sys.exit(1)

	return (filename, coil)



###############################################################################
# MAIN
###############################################################################

def main():
	usage = "%prog [options] --output-dir OUTPUTDIR MINCFILE [MINCFILE..]\n" + 9*" " + "-h for help"
	p = optparse.OptionParser(usage=usage, version="%prog " + version)
	p.add_option('--output-dir', '-o',                   dest="outdir",        default="discorr", help="output directory (default: discorr)")
	p.add_option('-n', '--dry-run', action="store_true", dest="dryrun",        default=False,     help="Print commands without executing them")
	p.add_option('--interp', '-i',                       dest="interp",        default="sinc",    help="resample interpolation (sinc,tricubic,trilinear)")
	p.add_option('--likefile', '-l',                     dest="likefile",      default={1:coil1_likefile,2:coil2_likefile,3:coil3_likefile}, help="likefiles")
	p.add_option('--queues', '-q',                       dest="queues",        default='',        help="queues (default: any)")
	p.add_option('--verbose', '-v', action="store_true", dest="verbose",       default=False,     help="Prints more information.")

	#p = argparse.ArgumentParser(description='Distortion correction for insert coil')
	#parser.add_argument('--output-dir', '-o', dest="outdir", default="discorr", help="output directory (default: discorr)")
	#parser.add_argument('--verbose', '-v', action="store_true", dest="verbose", default=False, help="Prints more information.")
	#parser.add_argument('mincfiles', metavar='file', type=string, nargs='+', help='MINC files to be distortion corrected')
	#args = parser.parse_args()

	options, arguments = p.parse_args()

	if len(arguments) < 1:
		p.error("no minc file to process")

	options.interp = '-' + options.interp

	for arg in arguments:
		(filename,coil) = splitarg(arg)
		if coil == coil and getcoil(filename)[1] != coil:
			print >> sys.stderr, 'WARNING: coil in file header "%s" does not match requested coil "%s"' % (getcoil(filename)[1],coil)
		
		print 'sge_batch', options.queues, '-l vf=4G -J dc-4G-rot36_b0 mincresample -2', options.interp, 'transform /projects/mice/matthijs/distortion-correction-november-2007/coil3-phantom-16oct2007/coil3-mri-to-CAD-lsq6-nlin-extra-large-super-inverted-lsq6.xfm -like ../raw/rot03_b0.mnc', filename, 'discorr/rot36_b0.november_2007_distortion_corrected.mnc'

if __name__ == '__main__':
	main()
	sys.exit(0)



#my $remove_registration = 0;
#
#      );
#GetOptions(\@arg_table, \@ARGV, \@left_over_args) or die "\n";
#
#my @mice = @left_over_args;
#
##if the -output-dir was not used, but the first argument
##is a directory, use this value to get the output directory
#if(!$output_dir and -d $mice[0])
#{
#	$output_dir = shift @mice;
#}
#
#if(!$output_dir)
#{
#	print "\n\nError: please specify an output directory.\n\n";
#	die $usage;
#}
#
#die $usage unless ($output_dir and $#mice >= 0);
#
#
#
## 1) Perform the correct transformation to the files
##
## The coil each file was scanned in is retrieved
## from the vnmr:coil attribute in the MINC header. 
#
#RegisterPrograms(["mincinfo", "mincresample", "sge_batch"]);
#system("mkdir -p $output_dir") unless (-d $output_dir);
#
#foreach my $mouse (@mice) 
#{
#	#check whether a coil number for this mouse was specified
#	my @parts = split(":", $mouse);
#	my $set_coil = undef;
#	if($parts[1] == 1 or $parts[1] == 2 or $parts[1] == 3)
#	{
#		$set_coil = $parts[1];
#		$mouse = $parts[0];
#	}
#	#if a coil was specified explicitly, but wasn't 1,2 or 3, exit:
#	if($#parts == 1 and !$set_coil)
#	{
#		print "\n\nError: the specified coil number should be 1,2 or 3. Was $parts[1]\n\n";
#		die;
#	}
#		
#	my ($dir, $base, $ext) = split_path($mouse, 'last', [qw(gz z Z)]);	
#	#at this point, $mouse should hold a minc file, verify this:
#	if(!($ext eq '.mnc'))
#	{
#		print "\n\nError: $mouse is not a minc file\n\n";
#		die;
#	}
#	my $output = $base;
#	$output .= ".november_2007_distortion_corrected.mnc";
#	$output = "${output_dir}/${output}";
#	
#	#transform options
## 	my @resample_opts = ("mincresample", $resampling, "-2", "-use_input_sampling");
#	my @resample_opts = ("mincresample", $resampling, "-2");
#	push @resample_opts, "-clobber", if $Clobber;
#	
#	my $coil;
#	
#	if(!$set_coil)
#	{
#		eval
#		{
#			Spawn(["mincinfo", "-attval", "vnmr:coil", $mouse], stdout=>\${coil});
#		};
#		if($@)
#		{
#			die "\nError: the vnmr:coil entry for file $mouse must be 1, 2, or 3 (was: $coil )\n\n";
#		}
#			$coil =~ s/.*(\d).*/$1/;
#	}
#	else
#	{
#		$coil = $set_coil
#	}
#	
#	if($coil == 1)
#	{
#		push @resample_opts, ("-transform", ${coil1_xfm}, "-like", ${coil1_likefile});
#	}
#	elsif($coil == 2)
#	{
#		push @resample_opts, ("-transform", ${coil2_xfm}, "-like", ${coil2_likefile});
#	}
#	elsif($coil == 3)
#	{
#		push @resample_opts, ("-transform", ${coil3_xfm}, "-like", ${coil3_likefile});
#	}
#	else
#	{
#		die "\nError: the vnmr:coil entry for file $mouse must be 1, 2, or 3 (was: $coil )\n\n";
#	}
#		
#	push @resample_opts, ($mouse, $output);
#	
#	if($Execute) 
#	{
#    if(!$sge or $spawn) 
#    {
#      Spawn([@resample_opts]);
#    }
#    else 
#    {
#		if(!$print)
#		{
#			Spawn(["sge_batch", "-q", "bigmem.q,all.q", "-l", "vf=4G", "-J", "dc-4G-$base", @resample_opts]);
#		}
#		else
#		{
#			print join(" ",@resample_opts), "\n"
#		}
#    }
# 	}
#  else 
#  {
#    print "@resample_opts \n";
#  }
#}
#
#exit;
#
