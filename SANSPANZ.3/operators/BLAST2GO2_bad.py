
from myoperator import BlockOperator
import sys

class BLAST2GO2(BlockOperator):
        """
        Implement blast2go
        """
        def __init__(self, glob,verbose=True):
                sys.stderr.write("# Init blast2go\n")
                self.glob=glob
		# two spreadsheets (input, output)
		[self.data,self.goclass_data]=glob.use_sheets(["data","goclass"])
		# define input column nicknames
		[self.qpid_col,self.status_col,self.golist_col,self.ecwlist_col]=self.data.use_columns(['qpid','status','GOclass','evidence_code_weight'])
		# sefine output column headers
		self.goclass_data.use_columns(["qpid","ontology","goid","desc", "BLAST2GO_1", "BLAST2GO_2"])
                # filter hit list; fetch go annotations and evidence-code weights
                [self.d,self.l]=glob.use_operators(['Filter','B2G'])
                # fetch annotations of go classes
		glob.use_online_dictionaries(["GOIDELIC"]) # to get ontology
		# method-specific parameter
		self.GO_weight=55


        def process(self,block):
                # sort SSRL
                if len(block)==0: return
                self.d.process(block)  # Filter status
		for row in block: self.l.process(row) # B2G is RowOperator; creates GOclass GOclass_count   evidence_code_weight columns

		# split GOclass and evidence_code_weight lists to entries in goclass table
		sum1={}
		sum2={}
		for row in block:
			if row[self.status_col]=='False': continue
			goidlist=row[self.golist_col].split(' ')
			ecwlist=row[self.ecwlist_col].split(' ')
			n=len(goidlist)
			for i in range(0,n):
				goid=goidlist[i]
				ecw=float(ecwlist[i])
				if not goid in sum1: 
					sum1[goid]=0
					sum2[goid]=0.0
				sum1[goid]+=1
				sum2[goid]+=ecw
		# write table with goid entries
		for goid in sum1.keys():
			# write sum of counts to BLAST20_1 and sum of wights ro BLAST2GO_2
			datarow=[row[self.qpid_col],self.glob.ontology[goid],goid,self.glob.godesc[goid],str(sum1[goid]),str(sum2[goid])]
			self.goclass_data.append_row(datarow)

