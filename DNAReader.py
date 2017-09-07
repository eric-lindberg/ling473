import time
import numpy as np
import sys

keyMap = {'G':0, 'g':0, 'A':1, 'a':1, 'T':2, 't':2, 'C':3, 'c':3}
noChildren = np.array([None, None, None, None]) # I hope that the junk yard a few blocks from here some day burns down

class TrieNode:
    """Manage the individual node in the Trie"""

    def __init__(self, label=None):
        self.label = label          # Key value used to find this node
        self.children = np.array([None, None, None, None])

        if self.label != None:
            self.label = self.label.upper()

    def __getitem__(self, item):
        return self.children[item]

    """Add a child node to the Trie"""
    def addChild(self, key):
        if key in keyMap:
            key_index = keyMap[key]
            node = TrieNode(key)
            self.children[key_index] = node
            return node
        return None

class Trie:
    """Top level controller for the Trie"""

    def __init__(self):
        self.head = TrieNode()

    # Add a new entry into the Trie
    def add(self, word):
        current_node = self.head        # Start at top of Trie
        word_finished = True

        # Look for existing matches in the Trie for this word
        for i in range(len(word)):
            if word[i] in keyMap and not np.array_equal(current_node.children, noChildren):
                keyIndex = keyMap[word[i]]
                node = current_node.children[keyIndex]
                if node != None :
                    current_node = current_node.children[word[i]]
                else:
                    word_finished = False
                    break
            else:
                word_finished = False
                break

        # From that point, create the new components of the word
        if not word_finished:
            while i < len(word):
                current_node = current_node.addChild(word[i])
                i += 1

    # Fine all matches to Trie values in the given text
    def find_matches(self, text):
        matches = list()        # Create a list of results as list of [offset, value]
        pn = self.head          # Start at the top of the Trie

        # Iterate over all letters in the text
        i = 0
        while i < len(text):
            if i > 0 and i % 100000 == 0:
                sys.stdout.write('.')
                sys.stdout.flush()
            letter = text[i]            # Get the next letter and it's lookup index
            if letter in keyMap and not np.array_equal(pn.children,noChildren): # Are there any child nodes?
                keyIndex = keyMap[letter]
                node = pn.children[keyIndex]
                if node != None:                            # Is there a child node for this letter?
                    pn = node
                    j = i+1
                    if np.array_equal(pn.children, noChildren): # Is that the last letter of a match?
                        matches.append([i, text[i:j]])
                    done_matching = False                   # DONE!

                    while not done_matching:                # On to the next letter in the possible match
                        if j >= len(text):
                            break
                        letter = text[j]

                        # Iterate over the remaining letters, looking for a full match
                        if letter in keyMap and not np.array_equal(pn.children, noChildren):
                            keyIndex = keyMap[letter]
                            node = pn.children[keyIndex]    # Get the next matching letter
                            if node != None:                # If there is a matching letter
                                pn = node
                                j += 1
                                if np.array_equal(pn.children, noChildren): # Matched the entire word
                                    matches.append([i, text[i:j]])
                                    done_matching = True
                            else:
                                done_matching = True        # Failed to match that word
                        else:
                            done_matching = True

            pn = self.head              # Done. Return to the top of the Trie
            i += 1

        return matches


def outputFileData(matches):
    for matchpair in matches:
        if len(matchpair[1]) > 0:
            print matchpair[0]
            for matchdata in matchpair[1]:
                print ('\t%07x' % matchdata[0]) + '\t' + matchdata[1].upper()

def outputMatchData(matches):
    matchdict = dict()
    for matchpair in matches:
        if len(matchpair[1]) > 0:
            filename = matchpair[0]
            for matchdata in matchpair[1]:
                keyWord = matchdata[1].upper()
                if keyWord not in matchdict:
                    matchdict[keyWord] = list()
                matchdict[keyWord].append([matchdata[0], filename])

    extra_credit_file = open('extra-credit', 'w+')
    for key in matchdict:
        extra_credit_file.write(key + '\n')
        for location_data in matchdict[key]:
            extra_credit_file.write(('\t%07x' % location_data[0]) + '\t' + location_data[1] + '\n')


def readTargetData():
    trie = Trie()
    # target_file = open('/opt/dropbox/17-18/473/project4/targets', 'r')
    target_file = open('/home/eric/Documents/ling473/project4/targets', 'r')
    for target_line in iter(target_file):
        trie.add(target_line.rstrip())

    return trie

def readAndMatchDataFiles(trie):

    file_indices = range(1,2)
    file_indices.append('X')
    file_indices.append('Y')
    total_start_time = time.time()
    all_matches = list()
    # directory_path = '/opt/dropbox/17-18/473/project4/hg19-GRCh37/'
    directory_path = '/home/eric/Documents/ling473/project4/'

    for file_index in file_indices:

        filename = 'chr{}.dna'.format(file_index)
        search_file = open(directory_path + filename, 'r')
        bigdata = search_file.read()

        start_time = time.time()
        all_matches.append([filename, trie.find_matches(bigdata)])

        end_time = time.time()

        print 'File ' + filename + ' completed in {} seconds'.format(end_time - start_time)

    print 'Total time elapsed is {} minutes'.format((end_time - total_start_time) / 60)

    outputFileData(all_matches)
    outputMatchData(all_matches)

# readAndMatchDataFiles(readTargetData())

def verifyOutput():
    trie = Trie()
    trie.add('CTGGAATATTCCCG')
    readAndMatchDataFiles(trie)

verifyOutput()