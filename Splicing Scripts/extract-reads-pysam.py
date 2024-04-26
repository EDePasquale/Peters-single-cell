#from __future__ import print_function
import pysam
import random
import sys
import os
import string
import fileinput

def getClusterBarcodes(barcodefile,reverseComplement=False):
	header=[]
	barcode_dict={}
	groups={}
	from Bio.Seq import Seq
	for line in open(barcodefile,'rU').xreadlines():
		line = line.rstrip()
		barcode,group_num,group_name=string.split(line,'\t')
		if '.' in barcode:
			barcode = string.split(barcode,'.')[0]
		#if '-1' in barcode:
		#	barcode = barcode[:-2]
		if reverseComplement:
			seq = Seq(barcode); barcode=seq.reverse_complement()
			barcode = str(barcode)
		if group_name in groups:
			output = groups[group_name]
		else:
			outfile = output_dir+'/'+ group_name+'.bam'
			output = pysam.AlignmentFile(outfile, "wb", template=bam)
			groups[group_name] = output
		barcode_dict[barcode]=output
	return barcode_dict,groups

def exportClusterBAMs(barcode_dict,groups):
	for read in bam.fetch():
		#print read.get_tag
		try: 
			barcode = read.get_tag('CB')
			#print barcode
			output = barcode_dict[barcode]
			output.write(read)
		except Exception: 
			pass

	bam.close()
	for group in groups:
		groups[group].close()
	
	output.close()

if __name__ == '__main__':
    import getopt
    reverse = False
    if len(sys.argv[1:])<=1:
		print 'Insufficient command-line arguments'; sys.exit()
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['i=','o=','f=','reverse='])
        #print sys.argv[1:]
        for opt, arg in options:
            if opt == '--i': input_file=arg
            elif opt == '--f': barcodefile = arg
            elif opt == '--reverse': reverse = True
            elif opt == '--o': output_dir = arg

	bam = pysam.AlignmentFile(input_file)
	barcode_dict,groups = getClusterBarcodes(barcodefile,reverseComplement=reverse)
	exportClusterBAMs(barcode_dict,groups)
    for file in os.listdir(output_dir):
        if '.bam' in file and '.bai' not in file:
            try: pysam.index(output_dir+'/'+file)
            except: pass
