"""Implementation of the Needleman-Wunsch algorithm for string alignment."""

import sys


_weights_file = 'SegmentDifferenceListForCatalanPolymorphemic.txt'

# dictionary of (phoneme, phoneme) pairs and their penalty scores
_penalty_scores = {}


def _load_scores ():
    """Calculate mismatch scores for each pair of phonemes according to the weights file."""

    with open(_weights_file, 'r') as wfile:
        next(wfile)

        # weights for each feature
        feature_weights = list(map(float, wfile.readline().split()))

        # set special gap penalty value
        _penalty_scores['-'] = feature_weights[0]

        for line in wfile:
            line = line.split()

            if line[0] == 'GAP' or line[1] == 'GAP':
                continue

            phonemes = tuple(line[0:2])
            feature_match = list(map(int, line[3:]))

            # calculate mismatch score
            _penalty_scores[phonemes] = sum([x * y for x, y in zip(feature_match, feature_weights)])


def _penalty (key):
    """
    Get the penalty score for an alignment of the two phonemes.

    Args:
        key (tuple): The phoneme pair

    Returns: The pre-calculated penalty score
    """

    if not _penalty_scores:
        _load_scores()

    try:
        return _penalty_scores[key]
    except KeyError:
        print(f'No penalty score could be calculated for the phoneme pair {key}.', file = sys.stderr)
        exit(1)


def align (lstr1, lstr2):
    """
    Calculate the optimal alignment of the strings according to the Needleman-Wunsch algorithm.

    Args:
        lstr1 (list): The first string, given as a list of phonemes
        lstr2 (list): The second string, given as a list of phonemes

    Returns: A tuple consisting of a tuple of the two aligned strings, given as a list of phonemes, and the calculated score
    """

    len1 = len(lstr1)
    len2 = len(lstr2)

    # initialize the first matrix row (all deletions)
    lastrow_scores = [_penalty('-') * i for i in range(len1 + 1)]
    lastrow_strings = [(lstr1[:i], ['-'] * i) for i in range(len1 + 1)]

    # calculate each subsequent matrix row
    for j in range(len2):
        # initialize the first row element (all insertions)
        row_scores = [_penalty('-') * (j + 1)]
        row_strings = [(['-'] * (j + 1), lstr2[:j + 1])]

        for i in range(len1):
            # calculate scores for matching, insertion, and deletion
            match_score = lastrow_scores[i] + _penalty((lstr1[i], lstr2[j]))
            ins_score = lastrow_scores[i + 1] + _penalty('-')
            del_score = row_scores[i] + _penalty('-')

            score = min(match_score, ins_score, del_score)
            row_scores.append(score)

            # derive new string alignment for the best scoring option
            if score == match_score:
                # matching: add next character for both strings
                row_strings.append((lastrow_strings[i][0] + [lstr1[i]], lastrow_strings[i][1] + [lstr2[j]]))
            elif score == ins_score:
                # insertion: add gap for string 1, next character for string 2
                row_strings.append((lastrow_strings[i + 1][0] + ['-'], lastrow_strings[i][1] + [lstr2[j]]))
            else:
                # deletion: add next character for string 1, gap for string 2
                row_strings.append((row_strings[i][0] + [lstr1[i]], row_strings[i][1] + ['-']))

        lastrow_scores = row_scores
        lastrow_strings = row_strings

    # optimal alignment: final element of the final row of the matrix
    return (row_strings[len1], round(row_scores[len1], 2))


if __name__ == '__main__':
    """Align two strings given as input."""

    if (len(sys.argv) != 3):
        print('Give two strings as arguments to find their optimal alignment.', file = sys.stderr)
        print('The phonemes of each should be separated with spaces.', file = sys.stderr)
        exit(1)

    strings, score = align(sys.argv[1].split(), sys.argv[2].split())

    # make phonemes line up for printing
    for i in range(len(strings[0])):
        phoneme_len = max(len(strings[0][i]), len(strings[1][i]))
        strings[0][i] = strings[0][i].ljust(phoneme_len)
        strings[1][i] = strings[1][i].ljust(phoneme_len)

    print(f'Alignment:')
    print('  ', ' '.join(strings[0]))
    print('  ', ' '.join(strings[1]))
    print()
    print(f'Score: {score}')
