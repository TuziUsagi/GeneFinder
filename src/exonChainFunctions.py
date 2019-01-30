import numpy as np
import pandas as pd
#this  function will find the highest score among the pathes (align, left gap, up gap),and assign the path to the pathMatrix 
# based on which path gives the highest score.
# the maximal score and path(s) are returned
def getMax(align, leftGap, upGap):
	scoreList = [align, leftGap, upGap]
	# each cell in pathMatrix stores a list of three 0/1 number, called path. Each number works as a boolean flag for path.
	# path[0] is the flag for path from upper left cell, path[1] flags path from left, path[2] flags the path from right. 
	# e.g. if the path is from the left and upper left, but not from the upper cell, the list path would be [1,1,0].
	maxScore = np.amax(scoreList, axis=0)
	path = [0,0,0]
	if maxScore <=0:
		maxScore = 0
		path = [0,0,0]
	else:
		for i in range(0,3):
			if scoreList[i] == maxScore:
				path[i] = 1
	return maxScore, path       

# find all the pathes that end at the location where the score is the maximal at the scoreMatrix
# a list of dictionary named allAlignment will be returned. Dictionary in allAlignment contains two keys: 'alignedSeq_1' and 'alignedSeq_2'
# the aligned sequences are constructed and stored in he allAlignment

# The search for all pathes is implemented in a recursive way.
# recursion for finding all the pathes that end at pathMatrix[rowIndex, colIndex]
def getAlignment(seq_1_split, seq_2_split, pathMatrix, rowIndex, colIndex, alignment):
	pathDirection = pathMatrix[rowIndex][colIndex]
	# if the search reaches the begining of a path, create a new record in the dictionary allAlignment
	if pathDirection == [0,0,0]:
		return alignment
	# recursive calls based on the type of path stored in the pathMatrix cell
	elif pathDirection == [1,0,0] or [1,1,0] or [1,1,1] or [1,0,1]:
		alignment = {'row_idx':rowIndex,'col_idx':colIndex}
		alignment = getAlignment(seq_1_split, seq_2_split, pathMatrix, rowIndex-1, colIndex-1, alignment)
		return alignment
	elif pathDirection == [0,1,0] or [0,1,1]:
		alignment = getAlignment(seq_1_split, seq_2_split, pathMatrix, rowIndex, colIndex-1, alignment)
		return alignment
	elif pathDirection == [0,0,1]:
		alignment = getAlignment(seq_1_split, seq_2_split, pathMatrix, rowIndex-1, colIndex, alignment)
		return alignment

# print out all the possible alignments, and save all alignment locally as 'alignment.txt'.
def generateOutput(scoreMatrix, pathMatrix,seq_1_split, seq_2_split, highestScore, highestIndex):    
	highestRow = highestIndex[0]
	highestCol = highestIndex[1]
	alignment = getAlignment(seq_1_split, seq_2_split, pathMatrix, highestRow, highestCol,[])
	alignment['highestCol']=highestCol         
	return alignment

def align_local (seq_1, seq_2, scoreParameters,genomeStart):
	dim_1 = len(seq_1) + 1 # the dimension for rows in the DP matrix
	dim_2 = len(seq_2) + 1 # the dimension for columns in the DP matrix
	# initiate scoreMatrix and pathMatrix
	scoreMatrix = [[0 for j in range(0,dim_2)] for i in range(0,dim_1)] # create score matrix of zeros with the defined dimensions from input sequences
	pathMatrix = [[[0,0,0] for j in range(0,dim_2)] for i in range(0,dim_1)] # create path matrix of zeros with the defined dimensions from input sequences
	# define the score for match, mismatch and gap penalty
	matchScore= scoreParameters[0]
	mismatchScore = scoreParameters[1]
	gap = scoreParameters[2]  
	# split the input sequences into a single characters
	seq_1_split = list(seq_1)
	seq_2_split = list(seq_2)
	# initiate the highestScore as 0
	highestScore = 0
	# build scoreMatrix and pathMatrix
	for row in range(1, dim_1):
		for col in range(1, dim_2):
			if seq_1_split[row-1] == seq_2_split[col-1]:
				align = scoreMatrix[row-1][col-1] + matchScore
			else:
				align = scoreMatrix[row-1][col-1] + mismatchScore
			leftGap = scoreMatrix[row][col-1] + gap
			upGap = scoreMatrix[row-1][col] + gap
			# find the maximal score among align, leftGap, and upGap, get the path with maximal score
			score, path = getMax(align, leftGap, upGap)
			scoreMatrix[row][col] = score
			pathMatrix[row][col] = path
			# replace the highestScore with current score if the current score is higher than the highestScore
			# store the location of the highestScore in highestIndex
			# if highestScore appeared multiple times, store the location for all of them
			if score >= highestScore:
				highestScore = score
				highestIndex = [row,col]
            
	# use the scoreMatrix, pathMatrix, and the location of the highest score to generate all the alignments 
	# seq_1_split, seq_2_split are passed for constructing the aligned sequences
	alignment = generateOutput(scoreMatrix, pathMatrix,seq_1_split, seq_2_split, highestScore, highestIndex)
	summary = {'genome_start': alignment['col_idx']+1+genomeStart, 'genome_end':alignment['highestCol']+genomeStart, 'score':highestScore}
	return summary

#calculate the length of a putative exons
def calLength(row):
	length = row['genome_end']-row['genome_start']
	return length

#calculate the rank of exon locations
def getRank(inputs):
	numOfIntervals = len(inputs['length'])
	intervalID = list(range(0,numOfIntervals)) + list(range(0,numOfIntervals))
	locs = inputs['genome_start']+inputs['genome_end']
	labels = ['start_rank']*numOfIntervals+['end_rank']*numOfIntervals
	score = inputs['length'] + inputs['length']
	allLocs = {'ID':intervalID,'loc': locs ,'label':labels, 'score': score}
	ranks = pd.DataFrame(data = allLocs)
	ranks['loc_Rank'] = ranks['loc'].rank(method='first')
	return ranks

# define class Stack
class Stack:
	def __init__(self):
		self.items = []
	def isEmpty(self):
		return self.items == []
	def push(self, item):
		self.items.append(item)
	def pop(self):
		return self.items.pop()
	def peek(self):
		return self.items[len(self.items)-1]          
