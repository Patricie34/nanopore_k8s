#!/bin/bash
# Download the CSV file from a URL
wget "https://docs.google.com/spreadsheets/d/1Rlglf6H8WmyEVRrKaOvDd5cHeYJaw-qqBEs8PMnaYQc/export?format=csv" -O "samplesheet.csv"

# Convert the CSV file to JSON using a Python script
python CsvToJson.py samplesheet.csv samplesheet.json
