#!/usr/bin/env python
import argparse
import yaml
import os
import subprocess
import shutil, errno
import sys

# we want to be agnostic to where the script is ran
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
WORKER_PATH = os.path.realpath(os.curdir)

# Function to copy the output directory content
def copyfolder(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise
        
# read the task definition yaml file
with open(os.path.join(SCRIPT_PATH, "rbpprofiles.yaml"), "r") as task_f:
	task_definition = yaml.load(task_f)

parser = argparse.ArgumentParser(
   description='Launches RBPprofiles with properly parset parameters')

# parser.add_argument(
#    '-fileA', type=str, default=["none"], nargs=1, help='Dataset A')
# 
# parser.add_argument(
#    '-fileB', type=str, default=["none"], nargs=1, help='Dataset B')

parser.add_argument(
   '-output_dir', type=str, nargs=1,
   help='Directory where the output is going to be stored')

# accept form fields
for item in task_definition['form_fields']:
   nargs = 1 if item['required'] else "?"
   parser.add_argument(
       '-FORM%s'%item['name'], type=str, default=['none'], nargs=nargs,
       help='Form argument: %s' % item)

# this parse_args stops if any unexpected arguments is passed
args = parser.parse_args()

OUTPUT_PATH = os.path.join(WORKER_PATH, args.output_dir[0])

# import code; code.interact(local=locals())
# import IPython
# IPython.embed()

import re
import StringIO
from Bio import SeqIO

Ppat = re.compile('>.*?\n[ARNDCQEGHILKMFPSTWYV]+', re.IGNORECASE)
if Ppat.match(args.FORMprotein_seq[0]) == None:
	args.FORMprotein_seq[0] = ">input_protein\n"+args.FORMprotein_seq[0]
protSeq = []
for record in SeqIO.parse(StringIO.StringIO(args.FORMprotein_seq[0]), "fasta"):
	protSeq.append(record)

protFile = os.path.join(OUTPUT_PATH,"protein.fasta")
output_handle = open(protFile, "w")
# SeqIO.write(protSeq, output_handle, "fasta")
output_handle.write(">input_protein\n"+record.seq)
output_handle.close()

#import IPython
#IPython.embed()

os.chdir(SCRIPT_PATH)
# print(WORKER_PATH)

random_number = (""" "{}" """.format(args.output_dir[0])).split("/")[3]

#import IPython
#IPython.embed()
if type(args.FORMtitle)==list:
	title = "none"
else:
	title = args.FORMtitle.replace(" ", "_")

# import IPython
# IPython.embed()

cmd = """bash runRBProfiles.sh "{}" "{}" "{}" "{}" """.format(random_number, args.FORMemail, title, protFile )

p = subprocess.Popen(cmd, cwd=SCRIPT_PATH, shell=True)

p.wait()

if p.returncode == 0:

	TMP_PATH = "./tmp/{}/outputs/".format(random_number)
	dirList=os.listdir(TMP_PATH)
	for file in dirList:
		shutil.copyfile(TMP_PATH+file, OUTPUT_PATH+file)
	if os.path.exists(OUTPUT_PATH+"images") == False :
		copyfolder(SCRIPT_PATH+"/images", OUTPUT_PATH+"images")
	
	from django.template import Template
	from django.template import Context
	from django.conf import settings
	from django.template import Template
	
	settings.configure(TEMPLATE_DIRS=(os.path.join(SCRIPT_PATH,'./')), DEBUG=True, TEMPLATE_DEBUG=True)

	with open(os.path.join(OUTPUT_PATH,"RBP.out"), "r") as rbpout:
		rbp_lines = rbpout.readlines()
		rbp_pfam  = rbp_lines[1].replace("#domain: ", "")
		rbp_scale = rbp_lines[2].replace("#scale:", "")
	
	# read the template file into a variable
	with open(os.path.join(SCRIPT_PATH, "index.html"), "r") as template_file:
	   template_string = "".join(template_file.readlines())
	
	import datetime
	
	# create template from the string
	t = Template(template_string)
	
	# context contains variables to be replaced
	c = Context(
	   {
		   "title": title,
		   "proteinFile" : protFile,
		   "RBPpfam" : rbp_pfam,
		   "RBPscale" : rbp_scale,
		   "randoms" : random_number,
		   "generated" : str(datetime.datetime.now()),
	   }
	)
	
	# and this bit outputs it all into index.html
	with open(os.path.join(OUTPUT_PATH, "index.html"), "w") as output: 
	   output.write(t.render(c))

else:
	sys.exit("The execution of the bash script failed.")
