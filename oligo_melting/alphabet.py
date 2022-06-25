"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from typing import Any, Dict, Iterable, Iterator, List, Tuple, Union


class AlphabetBase:
    __letters: Any

    def __init__(self, letters: str):
        if not letters:
            raise AssertionError("cannot build an empty alphabet")
        self.__letters = letters

    def __iter__(self) -> Iterator[str]:
        return iter(self.__letters)


class AlphabetBytes(AlphabetBase):
    __is_bytes: bool
    __data: Dict[bytes, str]

    def __init__(self, letters: str):
        super(AlphabetBytes, self).__init__(letters)
        self.__is_bytes = len(letters) <= 4
        if self.__is_bytes:
            self.__build_bytes()

    @property
    def is_bytes(self) -> bool:
        return self.__is_bytes

    def __build_bytes(self) -> None:
        pass
