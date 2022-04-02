import subprocess
import os
import sys
import shutil
# import progressbar


def input_file(file_name):
	with open(file_name, 'r') as file:
		file_content = file.readlines()
	return file_content


def output_file(file_name, content):
	with open(file_name, 'w') as file:
		file.writelines(content)


def RunCmd(cmd):
	pro = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# print(pro.stderr.read())
	pro.wait()


def MakeFtp(assembly):
	genomePair = {}
	for line in input_file(assembly):
		acc = line.strip().split('\t')[0]
		try:
			ftp = line.strip('\n').split('\t')[19] + '/' + line.strip('\n').split('\t')[19].split('/')[-1] + '_genomic.fna.gz'
			genomePair[acc] = ftp
		except IndexError:
			pass
	return genomePair


def ReadGenomeAccession(GenomeInfFile, assemblyPair, DownloadedDir, OutDir):
	fnaDir = os.path.join(OutDir, 'fna')
	try:
		os.mkdir(fnaDir)
	except FileExistsError:
		pass
	GenomesList = []
	for line in input_file(GenomeInfFile):
		genomeAcc = line.strip().split('\t')[2]
		if 'GCA' in genomeAcc:
			genomeAcc = genomeAcc
			genomeFile = os.path.join(DownloadedDir, '%s.fna' % genomeAcc)
			ToFile = os.path.join(fnaDir, '%s.fna' % genomeAcc)
			if os.path.isfile(genomeFile):
				shutil.copyfile(genomeFile, ToFile)
			else:
				ftpAdress = assemblyPair[genomeAcc]
				Cmd = ' '.join(['wget', '-c', ftpAdress, '-O', '%s.gz' % ToFile])
				RunCmd(Cmd)
				os.system('gzip -d %s.gz' % ToFile)
		else:
			genomeAcc = 'IMG%s' % genomeAcc
			genomeFile = os.path.join(DownloadedDir, '%s.fna' % genomeAcc)
			ToFile = os.path.join(fnaDir, '%s.fna' % genomeAcc)
			if os.path.isfile(genomeFile):
				shutil.copyfile(genomeFile, ToFile)
			else:
				print(genomeAcc)
		GenomesList.append(genomeAcc)
	return GenomesList


def RunSpecI(inputDir, Genome, threads, OutGenomeDir):
	inputFna = os.path.join(inputDir, Genome, '%s.ffn' % Genome)
	inputFaa = os.path.join(inputDir, Genome, '%s.faa' % Genome)
	outPut = os.path.join(OutGenomeDir, Genome)
	script = '/media/disk2/wangmin_Streptomyces/AM_synomys/SpecI/specI_v1.0/MOCATFetchMGs03.pl'
	cmd = ' '.join(
		['perl', script, '-v', '-o', outPut, '-t', threads, '-c', 'all', '-m', 'extraction', '-d', inputFna, inputFaa])
	RunCmd(cmd)


def GenomeAnnotation(InputDir, InputGenome, threads, ProkkaDir):
	InputFile = os.path.join(InputDir, '%s.fna' % InputGenome)
	OutDir = os.path.join(ProkkaDir, InputGenome)
	cmd = ' '.join(
		['prokka', "--outdir", OutDir, "--prefix", InputGenome, "--locustag", "'%s|ORF'" % InputGenome, "--noanno",
		 "--norrna", "--notrna", "--cpus", threads, InputFile])
	RunCmd(cmd)


def ExtractMarkerGenes(GenomesPath, GenomesList, OutputDir, RunThreads):
	ProkkaDir = os.path.join(OutputDir, 'prokka')
	MarkerGenesDir = os.path.join(OutputDir, 'MarkerGenes')
	HitTableDir = os.path.join(OutputDir, 'HitTable')
	FaaDir = os.path.join(OutputDir, 'proteins')
	try:
		os.mkdir(ProkkaDir)
	except FileExistsError:
		pass
	try:
		os.mkdir(HitTableDir)
	except FileExistsError:
		pass
	try:
		os.mkdir(FaaDir)
	except FileExistsError:
		pass
	try:
		os.mkdir(MarkerGenesDir)
	except FileExistsError:
		pass
	# progressbar_widgets_set = [
	# 	'Extract Marker Genes from Genomes', progressbar.Percentage(), progressbar.Bar('#'), progressbar.Timer(), ' | ',
	# 	progressbar.ETA()]
	# bar = progressbar.ProgressBar(widgets=progressbar_widgets_set, maxval=len(GenomesList))
	# bar.start()
	for num, genome in enumerate(GenomesList):
		genomeFile = os.path.join(GenomesPath, '%s.fna' % genome)
		# bar.update(num)
		if os.path.isfile(genomeFile):
			GenomeAnnotation(GenomesPath, genome, RunThreads, ProkkaDir)
			genomeFaa = os.path.join(ProkkaDir, genome, '%s.faa' % genome)
			RunSpecI(ProkkaDir, genome, RunThreads, MarkerGenesDir)
			genomeTable = os.path.join(MarkerGenesDir, genome, '%s.all.marker_genes_scores.table' % genome)
			shutil.copyfile(genomeFaa, os.path.join(FaaDir, '%s.faa' % genome))
			shutil.copyfile(genomeTable, os.path.join(HitTableDir, '%s.table' % genome))
	# bar.finish()


if __name__ == '__main__':
	GenomeInf = sys.argv[1:][0]
	assemblyFile = sys.argv[1:][1]
	genomeDownloaded = sys.argv[1:][2]
	runThreads = sys.argv[1:][3]
	outPut = sys.argv[1:][4]
	try:
		os.mkdir(outPut)
	except FileExistsError:
		pass
	assembly_pair = MakeFtp(assemblyFile)
	Genomes = ReadGenomeAccession(GenomeInf, assembly_pair, genomeDownloaded, outPut)
	GenomePath = os.path.join(outPut, 'fna')
	ExtractMarkerGenes(GenomePath, Genomes, outPut, runThreads)
