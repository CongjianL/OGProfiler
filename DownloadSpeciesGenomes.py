import subprocess
import os
import sys
import time


def input_file(file_name):
	with open(file_name, 'r') as file:
		file_content = file.readlines()
	return file_content


def output_file(file_name, content):
	with open(file_name, 'w') as file:
		file.writelines(content)


def split_line_inf(line):
	try:
		acc = line.strip('\n').split('\t')[0]
		genus = line.strip('\n').split('\t')[7].split(' ')[0].strip('[').strip(']')
		speices_name = line.strip('\n').split('\t')[7].split(' ')[1]
		Name = genus + ' ' + line.split('\t')[7].split(' ')[1]
		strains = line.strip('\n').split('\t')[8].split('strain=')[1]
		ftp = line.strip('\n').split('\t')[19] + '/' + line.split('\t')[19].split('/')[-1] + '_genomic.fna.gz'
		type_material = line.strip('\n').split('\t')[-1]
		release_time = line.strip('\n').split('\t')[14]
		assembly_statue = line.strip('\n').split('\t')[11]
		return Name, genus, speices_name, acc, strains, ftp, type_material, release_time, assembly_statue
	except IndexError:
		none_type = ''
		return none_type


def compare_time(time_forward, time_now):
	time1_f = time.strptime(time_forward, "%Y/%m/%d")
	time2_b = time.strptime(time_now, "%Y/%m/%d")
	if time1_f >= time2_b:
		return False
	else:
		return True


def search_assembly(assembly_file, targetSpecies):
	"""based on species name , we can search assemly list for getting genome information and fliting out it."""
	GenomeNum = 0
	assembly = input_file(assembly_file)
	speicesGenomeInf = {}
	for line in assembly[1:]:
		if split_line_inf(line) != '':
			Name, genus_name, speices_name, Acc, strains, ftp, type_material, release_time, assembly_statue = split_line_inf(
				line)
			if Name == targetSpecies:
				speicesGenomeInf[Acc] = [Name, strains, Acc, type_material, assembly_statue, release_time, ftp]
				GenomeNum += 1
	print('Searching Finished!\n %d Genomes Found' % GenomeNum)
	return speicesGenomeInf


def DownloadGenomes(GenomeInf, DownloadDir):
	for acc, inf in GenomeInf.items():
		FTP = inf[-1]
		path = os.path.join(DownloadDir, '%s.fna.gz' % acc)
		cmd = ' '.join(['wget', '-c', FTP, '-O', path, '-q'])
		process = subprocess.Popen(cmd, shell=True)
		process.wait()


if __name__ == '__main__':
	SearchSpecies = sys.argv[1:][1]
	assembly = sys.argv[1:][2]
	DownloadDir = sys.argv[1:][3]
	SeqInfName = sys.argv[1:][4]
	speicesGenomeInf_str = ''
	speicesGenomeInf = search_assembly(assembly, SearchSpecies)
	for name, inf in speicesGenomeInf.items():
		speicesGenomeInf_str += '\t'.join(inf[:-1]) + '\n'
	DownloadGenomes(speicesGenomeInf, DownloadDir)
	output_file(SeqInfName, speicesGenomeInf_str)
