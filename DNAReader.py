import time

class TrieNode:
    """Manage the individual node in the Trie"""
    def __init__(self, label=None):
        self.label = label          # Key value used to find this node
        self.children = dict()      # Next level entries in the Trie

        if self.label != None:
            self.label = self.label.upper()

    """Add a child node to the Trie"""
    def addChild(self, key):
        newNode = TrieNode(key)         # Create the new entry
        self.children[key] = newNode    # Add the new entry to the child list
        return newNode                  # Return this for functional call

class Trie:
    """Top level controller for the Trie"""
    def __init__(self):
        self.head = TrieNode()

    # Add a new word to the Trie
    def add(self, word):
        current_node = self.head        # Start at the top level of the Trie
        word_finished = True

        # First run through any existing matches for the first part of the word
        for i in range(len(word)):
            if word[i] in current_node.children:    # A match, so advance to the next level
                current_node = current_node.children[word[i]]
            else:
                word_finished = False               # No match, so start adding new letters
                break

        # Having advanced to the point where new letters are needed, add those
        if not word_finished:
            while i < len(word):    # Add letters for the remainder of the word
                current_node = current_node.addChild(word[i])
                i += 1

    # Find all matches to words in the Trie against the given text
    def find_matches(self, text):
        matches = list()        # Matches list, consisting of [offset, value] pairs

        pn = self.head          # Start at the top

        # Iterate over every letter in the text
        i = 0
        while i < len(text):
            letter = text[i].upper()
            if letter in pn.children :          # Is there a possible match?
                pn = pn.children[letter]        # Go to the next level and try to match
                j = i+1
                if len(pn.children) == 0:       # Adding check for corner case of one-character word
                    matches.append([i, text[i:j].upper()])
                done_matching = False
                while not done_matching:
                    if j >= len(text):          # Check for intervening end of file
                        break

                    letter = text[j].upper()

                    if letter in pn.children:   # If it still matches, keep on checking
                        pn = pn.children[letter]
                        j += 1
                        if len(pn.children) == 0:   # On full match, store the result
                            matches.append([i, text[i:j].upper()])
                            done_matching = True    # ... and go on to the next first letter
                    else:
                        done_matching = True        # Or if there's no match, we're also done

                pn = self.head                      # go back to the top of the Trie
            i += 1                  # Next first letter

        return matches


# Output data organized by files (output to stdout)
def outputFileData(matches):        # Data is sent in as a list of [offset, value] pairs
    for matchpair in matches:
        if len(matchpair[1]) > 0:
            print matchpair[0]
            for matchdata in matchpair[1]:
                print ('\t%07x' % matchdata[0]) + '\t' + matchdata[1]

# __EXTRA-CREDIT__
# Output data organized by value (output to extra-credit)
def outputMatchData(matches):       # Data is sent as a list of [offset, value] pairs
    matchdict = dict()
    for matchpair in matches:       # Iterate over the data and construct a dictionary of values
        if len(matchpair[1]) > 0:
            filename = matchpair[0]
            for matchdata in matchpair[1]:
                if matchdata[1] not in matchdict:
                    matchdict[matchdata[1]] = list()
                matchdict[matchdata[1]].append([matchdata[0], filename])

    # Output the resulting dictionary to the extra-credit file
    extra_credit_file = open('extra-credit', 'w+')
    for key in matchdict:
        extra_credit_file.write(key + '\n')     # The value that was matched
        for location_data in matchdict[key]:    # Output offset and filename pairs
            extra_credit_file.write(('\t%07x' % location_data[0]) + '\t' + location_data[1] + '\n')

# Manufacture the Trie of targets
def readTargetData():
    trie = Trie()
    target_file = open('/opt/dropbox/17-18/473/project4/targets', 'r')
    for target_line in iter(target_file):
        trie.add(target_line.rstrip())

    return trie

# Read one file and match it against the trie
def readAndMatchOneFile(trie, filename):
    directory_path = '/opt/dropbox/17-18/473/project4/hg19-GRCh37/'
    search_file = open(directory_path + filename, 'r')
    return [filename, trie.find_matches(search_file.read())]    # Returns the list of matches and offsets

# Read all chr*.dna files and run the contents against the trie
def readAndMatchDataFiles(trie):
    file_indices = range(1,23)      # chr1.dna, chr2.dna, ..., chr22.dna, chrX.dna, chrY.dna
    file_indices.append('X')
    file_indices.append('Y')
    all_matches = list()            # List of filename/matches pairs

    for file_index in file_indices:
        filename = 'chr{}.dna'.format(file_index)
        all_matches.append(readAndMatchOneFile(trie, filename))

    outputFileData(all_matches)     # Output the data by file
    outputMatchData(all_matches)    # Extra Credit - Output the matches by value to extra-credit file

readAndMatchDataFiles(readTargetData())

# DEBUG - Verifying matching algorithm against a single simple pattern on the real DNA data
# Also time the process to get a benchmark of time to complete
def verifyOutput():
    start_time = time.time()
    trie = Trie()
    trie.add('CTGGAATATTCCCG')
    readAndMatchDataFiles(trie)
    print 'Completed operation in {} minutes'.format((time.time() - start_time) / 60)

# verifyOutput()

# DEBUG - Verifying the Trie construction and matching algorithm on very small data
def testTrie():
    trie = Trie()
    trie.add('GATTACCA')
    trie.add('TAGACC')
    trie.add('CGTAA')

    print trie.find_matches('ATTAGATTACCATAGACCTAA')
    print trie.find_matches('NNNNNNNNNNNNNNNGATTACCAACGTAANNNNNNN')

# testTrie()

