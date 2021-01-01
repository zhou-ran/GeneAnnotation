#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/21 3:50 PM
__author__ = 'Zhou Ran'

"""
A interval object for deal with the genomic data.
discard for a while, there is no requirement to process the interval
"""
import math


def bs(arr, threshold):
    """
    Binary Search

    :param arr:
    :param threshold:
    :return:
    """
    low = 0
    high = len(arr)

    while low < high:
        mid = math.floor((low + high) / 2)

        if arr[mid] == threshold:
            return mid
        elif arr[mid] < threshold and mid != low:
            low = mid
        elif arr[mid] > threshold and mid != high:
            high = mid
        else:
            high = low = low + 1

    return low


class Interval:
    def __init__(self, lst):
        """
        User can input the [nested] tuple or [nested] list.
        You can add any custom information in lst[2:]
        :param lst:
        """
        self.interval = [i for i in Interval.__convert(lst)]

    @classmethod
    def __convert(cls, intervalinput):
        """
        convert the duty type into perfect type
        :param intervalinput:
        :return:
        """

        assert type(intervalinput) == list, "The input seem like not a list, {} found".format(type(intervalinput))

        if type(intervalinput) == list:
            return intervalinput

        elif type(intervalinput) == tuple:
            return map(lambda x: list(x), intervalinput)

        else:
            return [intervalinput]

    # @classmethod
    # def siteintersect(cls, site):



if __name__ == '__main__':
    a = Interval([])

    print(a.interval)

