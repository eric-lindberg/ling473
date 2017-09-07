class TrieNode:
    """Manage the individual node in the Trie"""
    def __init__(self, label=None):
        self.label = label          # Key value used to find this node
        self.children = dict()

        if self.label != None:
            self.label = self.label.upper()

    """Add a child node to the Trie"""
    def addChild(self, key):
        newNode = TrieNode(key)
        self.children[key] = newNode
        return newNode

class Trie:
    """Top level controller for the Trie"""

    def __init__(self):
        self.head = TrieNode()

    def add(self, word):
        current_node = self.head
        word_finished = True

        for i in range(len(word)):
            if word[i] in current_node.children:
                current_node = current_node.children[word[i]]
            else:
                word_finished = False
                break

        if not word_finished:
            while i < len(word):
                current_node.addChild(word[i])
                current_node = current_node.addChild(word[i])
                i += 1

    def find_matches(self, text):
        matches = list()

        pn = self.head

        for i in range(0, len(text)):
            letter = text[i].upper()
            if letter in pn.children :
                pn = pn.children[letter]
                j = i+1
                if len(pn.children) == 0:
                    matches.append([i, text[i:j].upper()])
                done_matching = False
                while not done_matching:
                    if j >= len(text):
                        break

                    letter = text[j].upper()

                    if letter in pn.children:
                        pn = pn.children[letter]
                        j += 1
                        if len(pn.children) == 0:
                            matches.append([i, text[i:j].upper()])
                            done_matching = True
                    else:
                        done_matching = True

                pn = self.head

        return matches


def outputFileData(matches):
    for matchpair in matches:
        if len(matchpair[1]) > 0:
            print matchpair[0]
            for matchdata in matchpair[1]:
                print ('\t%07x' % matchdata[0]) + '\t' + matchdata[1]

def outputMatchData(matches):
    matchdict = dict()
    for matchpair in matches:
        if len(matchpair[1]) > 0:
            filename = matchpair[0]
            for matchdata in matchpair[1]:
                if matchdata[1] not in matchdict:
                    matchdict[matchdata[1]] = list()
                matchdict[matchdata[1]].append([matchdata[0], filename])

    extra_credit_file = open('extra-credit', 'w+')
    for key in matchdict:
        extra_credit_file.write(key + '\n')
        for location_data in matchdict[key]:
            extra_credit_file.write(('\t%07x' % location_data[0]) + '\t' + location_data[1] + '\n')


def readTargetData():
    trie = Trie()
    target_file = open('/opt/dropbox/17-18/473/project4/targets', 'r')
    for target_line in iter(target_file):
        trie.add(target_line.rstrip())

    return trie

def readAndMatchDataFiles(trie):

    file_indices = range(1,23)
    file_indices.append('X')
    file_indices.append('Y')
    all_matches = list()
    directory_path = '/opt/dropbox/17-18/473/project4/hg19-GRCh37/'

    for file_index in file_indices:

        filename = 'chr{}.dna'.format(file_index)
        search_file = open(directory_path + filename, 'r')
        bigdata = search_file.read()

        all_matches.append([filename, trie.find_matches(bigdata)])

    outputFileData(all_matches)
    outputMatchData(all_matches)

readAndMatchDataFiles(readTargetData())

def verifyOutput():
    trie = Trie()
    trie.add('CTGGAATATTCCCG')
    readAndMatchDataFiles(trie)

# verifyOutput()

def testTrie():
    trie = Trie()
    trie.add('GATTACCA')
    trie.add('TAGACC')
    trie.add('CGTAA')

    print trie.find_matches('ATTAGATTACCATAGACCTAA')
    print trie.find_matches('NNNNNNNNNNNNNNNGATTACCAACGTAANNNNNNN')

# testTrie()
