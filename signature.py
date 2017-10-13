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

#pylog=open("./outputs/signaturepy.log","w")

# Function to copy the output directory content
def copyfolder(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise

if sys.argv[1] == "-text=Yes":
	# read the task definition yaml file
    with open(os.path.join(SCRIPT_PATH, "signature.yaml"), "r") as task_f:
        task_definition = yaml.load(task_f)
    input_mode = "text"
else:
    with open(os.path.join(SCRIPT_PATH, "signature_file.yaml"), "r") as task_f:
        task_definition = yaml.load(task_f)
    input_mode = "file"

parser = argparse.ArgumentParser(
   description='Launches catRAPID signature with properly parset parameters')

# parser.add_argument(
#    '-fileA', type=str, default=["none"], nargs=1, help='Dataset A')
#
# parser.add_argument(
#    '-fileB', type=str, default=["none"], nargs=1, help='Dataset B')
if input_mode == "file":
    parser.add_argument('-fileA', type=str, default=["none"], nargs=1, help='Fasta sequence')


parser.add_argument(
   '-output_dir', type=str, nargs=1,
   help='Directory where the output is going to be stored')
parser.add_argument(
    '-text', type=str, nargs="?", help='Just to enable text mode')

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


if input_mode == "text":
    Ppat = re.compile('>.*?\n[ARNDCQEGHILKMFPSTWYV]+', re.IGNORECASE)
    if Ppat.match(args.FORMprotein_seq[0]) == None:
    	args.FORMprotein_seq[0] = ">input_protein\n"+args.FORMprotein_seq[0]
    protSeq = []
    protFile = os.path.join(OUTPUT_PATH,"protein.fasta")
    output_handle = open(protFile, "w")
    for record in SeqIO.parse(StringIO.StringIO(args.FORMprotein_seq[0]), "fasta"):
    	protSeq.append(record)
    	output_handle.write(str(">input_protein\n"+record.seq))
    # SeqIO.write(protSeq, output_handle, "fasta")
    output_handle.close()

else:
    protSeq = []
    protFile = os.path.join(OUTPUT_PATH.replace("output/", ""),"protein.fasta")
    input_handle = open(args.fileA[0], "rU")
    for record in SeqIO.parse(input_handle, "fasta"):
        record.id=record.id.replace(" ","").replace("\t","")+record.description.replace(" ","").replace("\t","")
        record.description=""
        protSeq.append(record)
    output_handle = open(protFile, "w")
    SeqIO.write(protSeq, output_handle, "fasta")
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

cmd = """bash runsignature.sh "{}" "{}" "{}" "{}" """.format(random_number, args.FORMemail, title, protFile )

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

	import datetime

	#with open(os.path.join(OUTPUT_PATH,"noPass.txt"), "r") as rbpout:
#		rbp_lines = rbpout.readlines()
#		rbp_pfam  = rbp_lines[0].split()[0]
	#
	## read the template file into a variable

	try:
		#if jpeg is available
		test= open(os.path.join(OUTPUT_PATH,"plot.jpeg"), "r")
		#pylog.write("binding domains detected... ")
		with open(os.path.join(OUTPUT_PATH,"prediction.txt"), "r") as rbpout:
			rbp_lines = rbpout.readlines()

		PredictionScore=rbp_lines[1].split()[1]
                descr=rbp_lines[0][1:]

		c_score=rbp_lines[2].split()[1]
		nc_score=rbp_lines[3].split()[1]
		p_score=rbp_lines[4].split()[1]

		with open(os.path.join(OUTPUT_PATH,"HMMER.results"), "r") as hmmercheck:
			hmmer_lines=hmmercheck.readlines()


		if (float(PredictionScore)<0.5) and (descr=="onlyHMMER\n"):
			with open(os.path.join(SCRIPT_PATH, "indexonlyHMMERdetected.html"), "r") as template_file:
				   template_string = "".join(template_file.readlines())
			c = Context({"title": title, "proteinFile" : protFile, "PredictionScore" : PredictionScore, "classical": c_score, "Nonclassical" :nc_score, "putative":p_score, "randoms" : random_number, "generated" : str(datetime.datetime.now()), })

		elif (float(PredictionScore)>=0.5) and (hmmer_lines[0].split()[0]=="Nodomain"):
			with open(os.path.join(SCRIPT_PATH, "indexGammaNoTable.html"), "r") as template_file:
				   template_string = "".join(template_file.readlines())
			# context contains variables to be replaced
			c = Context({"title": title, "proteinFile" : protFile, "PredictionScore" : PredictionScore, "Description": descr, "classical": c_score, "Nonclassical" :nc_score, "putative":p_score, "randoms" : random_number, "generated" : str(datetime.datetime.now()), })
		else:
			with open(os.path.join(SCRIPT_PATH, "indexGammaWithTable.html"), "r") as template_file:
				   template_string = "".join(template_file.readlines())
			# context contains variables to be replaced
			c = Context({"title": title, "proteinFile" : protFile, "PredictionScore" : PredictionScore, "Description": descr, "classical": c_score, "Nonclassical" :nc_score, "putative":p_score, "randoms" : random_number, "generated" : str(datetime.datetime.now()), })




	except:

		##check if multiple submission is given
		#pylog.write("check if multiple submission... ")
		try:
			test= open(os.path.join(OUTPUT_PATH,"catRAPID_s_predictions.txt"), "r")
			#pylog.write("multiple submission... OK ")
			with open(os.path.join(SCRIPT_PATH, "multiSubmission.html"), "r") as template_file:
				template_string = "".join(template_file.readlines())

			c = Context({"title": title,"proteinFile" : protFile, "randoms" : random_number, "generated" : str(datetime.datetime.now()),})


		except:
		###else give error
			with open(os.path.join(SCRIPT_PATH, "indexError.html"), "r") as template_file:
				template_string = "".join(template_file.readlines())

			with open(os.path.join(OUTPUT_PATH,"prediction.txt"), "r") as rbpout:
				rbp_lines = rbpout.readlines()
				#rbp_pfam  = rbp_lines[0].split()[0]
				PredictionScore=rbp_lines[1].split()[1]
				descr=rbp_lines[0][1:-1]

			c = Context({"title": title,"proteinFile" : protFile, "PredictionScore" : PredictionScore,"randoms" : random_number,"Description": descr, "generated" : str(datetime.datetime.now()),})

		#c = Context({"title": title,"proteinFile" : protFile,"randoms" : random_number,"generated" : str(datetime.datetime.now()),})



	t = Template(template_string)
#	import datetime

	# create template from the string


	# context contains variables to be replaced
	#c = Context(
	 #  {
	#	   "title": title,
	#	   "proteinFile" : protFile,
	#	#   "RBPpfam" : rbp_pfam,
	#	 #  "RBPscale" : rbp_scale,
	#	   "randoms" : random_number,
	#	   "generated" : str(datetime.datetime.now()),
	 #  }
	#)

	# and this bit outputs it all into index.html
	with open(os.path.join(OUTPUT_PATH, "index.html"), "w") as output:
	   output.write(t.render(c))

else:
	sys.exit("The execution of the bash script failed.")
