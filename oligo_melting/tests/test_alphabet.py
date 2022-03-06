"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from oligo_melting.alphabet import AlphabetBase, AlphabetBytes


def test_AlphabetBase() -> None:
    alphabet_string = "ACTG"
    alphabet = AlphabetBase(alphabet_string)
    for letter_id, letter in enumerate(alphabet):
        assert alphabet_string[letter_id] == letter


def test_AlphabetBytes() -> None:
    alphabet_string = "ACTG"
    alphabet = AlphabetBytes(alphabet_string)
    assert alphabet.is_bytes
    for letter_id, letter in enumerate(alphabet):
        assert alphabet_string[letter_id] == letter
