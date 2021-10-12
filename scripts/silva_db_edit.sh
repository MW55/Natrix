#!/bin/bash
awk -F "\t" '{print $1"."$2"."$3"\t"$4}'  $1
