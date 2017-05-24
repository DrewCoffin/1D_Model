#!/bin/bash

sed -i "s/LNG_GRID\=[0-9]*/LNG_GRID=$1/g" varTypes.f90
