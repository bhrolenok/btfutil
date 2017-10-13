import glob, os.path, string, numpy, time
import tarfile
import numpy

VERBOSE_TIMEOUT=30 #Time in seconds between print outs

def verbose_readlines(infile):
	theline = infile.readline()
	lasttime = time.time()
	lastidx = 0
	curidx = 0
	while theline != "":
		yield theline
		theline = infile.readline()
		curtime = time.time()
		if curtime - lasttime > VERBOSE_TIMEOUT:
			if lastidx == 0:
				print "[BTFUtil] Loading",infile.name
			print "[BTFUtil] Line", curidx, "lps",(curidx-lastidx)/float(VERBOSE_TIMEOUT)
			lastidx = curidx
			lasttime = curtime
		curidx += 1

class BTF:
	def __init__(self,src_path=None):
		self.column_filenames = dict()
		self.column_data = dict()
		self.mask = None
		self.tfile = None
		if not(src_path is None):
			if os.path.isdir(src_path):
				self.import_from_dir(src_path)
			else:
				self.import_from_tar(src_path)

	def import_from_dir(self,dirname):
		self.mask = None
		new_columns = glob.glob(os.path.join(dirname,'*.btf'))
		for column in new_columns:
			cname = os.path.basename(column)[:-4]
			self.column_filenames[cname] = column
		self.tfile = None

	def import_from_tar(self,tarname):
		self.tfile = tarfile.open(tarname)
		new_columns = filter(lambda nme: nme[-4:].lower()=='.btf', self.tfile.getnames())
		for column in new_columns:
			cname = os.path.basename(column)[:-4]
			self.column_filenames[cname] = column

	def save_to_dir(self,dirname,overwrite=False):
		for key in self.column_data:
			fname = os.path.join(dirname,key+".btf")
			if os.path.exists(fname) and not(overwrite):
				raise IOError("File exists: ["+fname+"]")
			outf = open(fname,"w")
			for line in self.column_data[key]:
				outf.write(line+"\n")
			outf.close()

	def load_column(self,cname):
		if cname in self.column_filenames:
			if self.tfile is None:
				sourcef = open(self.column_filenames[cname])
			else:
				sourcef = self.tfile.extractfile(self.column_filenames[cname])
			self.column_data[cname] = tuple(map(string.strip, verbose_readlines(sourcef)))
			return True
		return False

	def load_columns(self,cnames):
		return {k: self.load_column(k) for k in cnames}

	def load_all_columns(self):
		return self.load_columns(self.column_filenames.keys())

	def has_columns(self,cnames):
		for c in cnames:
			if not(c in self):
				return False
		return True

	def filter_by_col(self,col,val=None):
		if val is None:
			self.mask = tuple(ele.capitalize() == 'True' for ele in self[col])
		else:
			self.mask = tuple(ele == val for ele in self[col])

	def to_nparr(self):
		for key in self.column_data:
			col_name = key
			col_data = numpy.array(map(lambda s: s.split(), self[key]))
			print col_data
			col_width = col_data.shape[1]
		rv = numpy.column_stack([self.column_data[key] for key in self.column_data])
		rv = numpy.array(rv,dtype=[(key,'float') for key in self.column_data])
		return rv

	def __contains__(self,key):
		return key in self.column_filenames

	def __getitem__(self,key):
		if not(self.__contains__(key)):
			raise KeyError("No column named ["+str(key)+"]")
		rv = None
		if key in self.column_data:
			rv = self.column_data[key]
		else:
			if self.load_column(key):
				rv = self.column_data[key]
			else:
				raise KeyError("Could not load BTF column: ("+str(key)+","+str(self.column_filenames[key])+")")
		if self.mask is None:
			return rv
		else:
			return filter(lambda d: not(d is None), tuple(map(lambda data,mask: data if mask else None, rv, self.mask)))

def timeseries(btf,fun,pColNames,tCol='clocktime'):
	oldT = None
	firstNewT = 0
	alldata = list()
	alltimes = list()
	num_items = len(btf[tCol])
	last_atime = time.time()
	start_atime = last_atime
	last_line = 0
	for tIdx in range(num_items):
		if (oldT is None):
			oldT = btf[tCol][tIdx]
		elif btf[tCol][tIdx] != oldT:
			cdata = numpy.column_stack([map(float,btf[pName][firstNewT:tIdx]) for pName in pColNames])
			alldata.append(fun(cdata))
			alltimes.append(float(oldT))
			firstNewT = tIdx
			oldT = btf[tCol][tIdx]
		cur_atime = time.time()
		if (cur_atime-last_atime)>VERBOSE_TIMEOUT:
			if last_atime == start_atime:
				print "[BTFUtil] timeseries"
			print "[BTFUtil]","%f%%"%(100.0*float(tIdx)/float(num_items)),"@",float(tIdx-last_line)/float(cur_atime-last_atime),"lps"
			last_atime=cur_atime
			last_line = tIdx
	return numpy.array(alltimes), numpy.array(alldata)

def printif(s,q):
	if q:
		print s

def split_subsequences(btf,subseq_length_t,ignore_shorter=True,depth=0,debug=False,frameBoundaryColName='timestamp'):
	done = False
	rv = tuple()
	lasttime = time.time()
	last_remaininglines = len(btf[frameBoundaryColName])
	total_lines = last_remaininglines
	while not(done):
		printif("depth %d"%depth,debug)
		head_btf = BTF()
		tail_btf = BTF()
		head_btf.column_filenames = btf.column_filenames
		tail_btf.column_filenames = btf.column_filenames
		seq_start_t = float(btf['clocktime'][0])
		block_start_idx = 0
		id_set = None
		max_len = len(btf[frameBoundaryColName])
		if ((float(btf['clocktime'][max_len-1])-seq_start_t)<subseq_length_t) and (ignore_shorter):
			printif("Final segment too short",debug)
			done = True
			break
		while block_start_idx < max_len and (float(btf['clocktime'][block_start_idx])-seq_start_t)<subseq_length_t:
			block_end_idx=block_start_idx
			tmp_id_set = set()
			while block_end_idx < max_len and float(btf[frameBoundaryColName][block_end_idx])==float(btf[frameBoundaryColName][block_start_idx]):
				tmp_id_set.add(btf['id'][block_end_idx])
				block_end_idx += 1
			if id_set is None:
				id_set = tmp_id_set
			block_start_idx = block_end_idx
			if id_set != tmp_id_set:
				break
		last_seq_idx = min(max_len-1,block_start_idx)
		printif("seq length %fs"%(float(btf['clocktime'][last_seq_idx])-seq_start_t),debug)
		printif("last_seq_idx %d"%last_seq_idx,debug)
		for key in btf.column_filenames:
			if not(key in btf.column_data):
				btf.load_column(key)
			tail_btf.column_data[key] = btf[key][last_seq_idx:]
			head_btf.column_data[key] = btf[key][:last_seq_idx]
		if (float(btf['clocktime'][last_seq_idx])-seq_start_t)<subseq_length_t and ignore_shorter:
			printif("ended early",debug)
			thing_to_add = tuple()
		else:
			thing_to_add = (head_btf,)
		rv += thing_to_add
		btf = tail_btf
		depth=depth+1
		curtime = time.time()
		if curtime - lasttime > VERBOSE_TIMEOUT:
			print "Remaining lines:",max_len,"({}%)".format(100.0*float(total_lines-max_len)/float(total_lines)),"lps:",float(last_remaininglines-max_len)/float(curtime - lasttime)
			lasttime = curtime
			last_remaininglines = max_len
	return rv

def writeInitialPlacement(outf,initialPlacementBTF,frameBoundaryColName='timestamp'):
	rowIdx = 0
	while rowIdx < len(initialPlacementBTF['id']) and initialPlacementBTF[frameBoundaryColName][rowIdx] == initialPlacementBTF[frameBoundaryColName][0]:
		outf.write(initialPlacementBTF['id'][rowIdx])
		outf.write(" "+initialPlacementBTF['xpos'][rowIdx])
		outf.write(" "+initialPlacementBTF['ypos'][rowIdx])
		outf.write(" "+initialPlacementBTF['timage'][rowIdx]+"\n")
		rowIdx += 1
	return rowIdx

def btf2data(btf,feature_names,augment,ys_colname='dvel'):
	features = numpy.column_stack([map(lambda line: map(float,line.split()), btf[col_name]) for col_name in feature_names])
	if augment:
		features = numpy.column_stack([features,numpy.ones(features.shape[0])])
	ys = numpy.array(map(lambda line: map(float, line.split()), btf[ys_colname]))
	return features,ys

def split_btf_trajectory(btf,feature_names,augment,id_colname='id'):
	features,ys = btf2data(btf,feature_names,augment)
	npid = numpy.array(map(int,btf[id_colname]))
	unique_ids = set(npid)
	return {eyed:features[npid==eyed] for eyed in unique_ids}

def merge_by_column(btf1,btf2,colname):
	merged = BTF()
	merged.column_filenames = {cname:None for cname in btf1.column_filenames.keys()}
	combo = btf1[colname]+btf2[colname]
	sorted_indexes = tuple(idx for idx,key in sorted(enumerate(combo),key=lambda x:x[1]))
	for cname in merged.column_filenames.keys():
		combo = btf1[cname]+btf2[cname]
		merged.column_data[cname] = tuple(combo[idx] for idx in sorted_indexes)
	return merged

def load_sequence_dir(seqdir):
	return [BTF(os.path.join(seqdir,name)) for name in os.listdir(seqdir) if os.path.isdir(os.path.join(seqdir,name))]
