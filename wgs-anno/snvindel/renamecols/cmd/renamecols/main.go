package main

import (
	"fmt"
	"os"
	"renamecols/internal/renameprocess"
	"strings"
)

var bedColName string = "HGMD_AllMut_202303"
var bed2ColName string = "HGMD_AllMut-collapse"

func main() {
	// Check number of args
	if len(os.Args) < 2 {
		fmt.Println("Please enter the input file path.")
		os.Exit(1)
	}

	// #1. Get sample names from VCF
	vcfPath := strings.TrimSuffix(os.Args[1], ".txt") + ".vcf"
	var samples []string = renameprocess.GetSamples(vcfPath)

	// #2. Rename columns in txt
	// Rename columns and return lines to write a new txt
	txtPath := os.Args[1]
	var lines []string = renameprocess.ProcessTxtFile(txtPath, samples, bedColName, bed2ColName)

	// #3. Write a new txt
	outPath := strings.TrimSuffix(txtPath, ".txt") + ".renamed.txt"
	renameprocess.SaveRenamedFile(lines, outPath)
}
