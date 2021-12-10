# transform target map according to the rotation and translation of each of the top 10 models in VESPER output
import os
import argparse
import sys
import Bio.PDB as bpdb
import numpy as np

def write_chimera_command(command_file, map1, map2, r_vector, t_vector, output_name):
	command_file = open(command_file, 'w')
	command_file.write('from chimera import runCommand as rc\n\n')
	command_file.write('rc("open ' + map1 + '")\n')
	command_file.write('rc("open ' + map2 + '")\n')

	# rotate map
	command_file.write('rc("turn z ' + r_vector[0] + ' center 0,0,0 coord #1 models #1")\n')
	command_file.write('rc("turn y ' + r_vector[1] + ' center 0,0,0 coord #1 models #1")\n')
	command_file.write('rc("turn x ' + r_vector[2] + ' center 0,0,0 coord #1 models #1")\n')

# 	command_file.write('rc("turn z ' + r_vector[0] + ' center #1 coord #1 models #1")\n')
# 	command_file.write('rc("turn y ' + r_vector[1] + ' center #1 coord #1 models #1")\n')
# 	command_file.write('rc("turn x ' + r_vector[2] + ' center #1 coord #1 models #1")\n')

	# translate map
# 	command_file.write('rc("move ' + ','.join(t_vector) + ' coord #1 model #1")\n')

	# save transformed map
# 	command_file.write('rc("vop #1 resample onGrid #0")\n')
# 	command_file.write('rc("volume #2 save ' + output_name + '")\n')
# 	command_file.write('rc("write relative #0 #1 ' + output_name + '")\n')
	command_file.write('rc("write #1 ' + output_name + '")\n')
	command_file.close()

# parse the arguments from command line
parser = argparse.ArgumentParser(description = 'Transform the target map given rotation and translation information in the VESPER output file.')
# parser.add_argument('-i1', '--input1', required = True, action = 'store', dest = 'ref_map', help = 'Required. Name of the reference map file.')
parser.add_argument('-i2', '--input2', required = True, action = 'store', dest = 'target_map', help = 'Required. Name of the target map file.')
parser.add_argument('-t', required = True, action = 'store', dest = 'vesper_result', help = 'Required. Name of the result file from VESPER.')
parser.add_argument('-output', required = True, action = 'store', dest = 'output_prefix', help = 'Prefix of the output files.')
parser.add_argument('-odir', action = 'store', dest = 'out_dir', help = 'Optional. Directory for the transformed target map files. If not specified, the transformed target map files would be written to the current directory')

args = parser.parse_args()

# ref_map = args.ref_map
target_map = args.target_map
chains = target_map.split(',')
print(chains)
# target = [float(i) for i in target_map.split(',')]
vesper_result = args.vesper_result
output_prefix = args.output_prefix

if not args.out_dir:
	out_dir = './'
else:
	out_dir = args.out_dir

command_file = 'chimera_command.py'

# extract rotation and translation information from vesper_result
for i in chains:
	transform_file = i+ vesper_result
	chain_file = out_dir + i+'.pdb'
	
	with open(transform_file) as result:
		#model_line_start = ['#0', '#1', '#2', '#3', '#4', '#5', '#6', '#7', '#8', '#9']
		model_num = 0
		for line in result:
			#if line[0:2] in model_line_start:
			if line[0].isdigit()==False:
				continue
			model_num += 1
			r_info = line.split()[1:4]
			t_info = line.split()[4:7]

			rotation_vector= [r_info[0], r_info[1], r_info[2]]
			translation_vector = [t_info[0], t_info[1], t_info[2]]
		
			print('rotation_vector ', rotation_vector)
			print('translation_vector ', translation_vector)
		
			#transform target map
			out_name = i + '_' + output_prefix + '_' + str(model_num) + '.pdb'
			write_chimera_command(command_file, chain_file, chain_file, rotation_vector, translation_vector, out_name)
			run_command = 'chimera --silent --nogui ' + command_file 
			os.system(run_command)

			#remove chimera command file
			os.system('rm ' + command_file)
			os.system('rm ' + command_file + 'c')
		
			pdbparser = bpdb.PDBParser(PERMISSIVE=1,QUIET = True)
			pdb_file = pdbparser.get_structure('',out_name)[0]
			center = np.average([atom.coord for atom in pdb_file.get_atoms()],axis=0)
			translation = [float(i) for i in translation_vector]

			translation -= center
			print('translation ', translation)

			ident_rot = np.eye(3)
			pdb_file.transform(ident_rot,(translation.T))

			io = bpdb.PDBIO()
			io.set_structure(pdb_file)
			io.save(out_name)
			#print('Result written to :',output_file)
