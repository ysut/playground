package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"

	"github.com/BurntSushi/toml"
)

type Config struct {
	RenameMap map[string]string `toml:"columns"`
}

func main() {
	// Check number of args
	if len(os.Args) < 3 {
		fmt.Println("Please enter an input file and a config file path.")
		os.Exit(1)
	}
	configPath := os.Args[2]

	var config Config
	if _, err := toml.DecodeFile(configPath, &config); err != nil {
		fmt.Fprintf(os.Stderr, "Error reading config file: %v\n", err)
		os.Exit(1)
	}

	// #1. Get sample names from VCF
	vcfPath := strings.TrimSuffix(os.Args[1], ".txt") + ".vcf"
	var samples []string = GetSamples(vcfPath)

	// #2. Rename columns in txt
	// Rename columns and return lines to write a new txt
	txtPath := os.Args[1]
	var lines = RenameProcess(txtPath, samples, config.RenameMap)

	// #3. Write a new txt
	outPath := strings.TrimSuffix(txtPath, ".txt") + ".renamed.txt"
	SaveRenamedFile(lines, outPath)
}

// Functions
func ChangeColName(colNames []string, renameMap map[string]string) []string {
	for i, colName := range colNames {
		for key, value := range renameMap {
			if colName == key {
				// Change column name to value
				colNames[i] = value
			}
		}
	}
	return colNames
}

func RenameProcess(txtPath string, samples []string, renamemap map[string]string) []string {
	// Open a file annotated by ANNOVAR
	txtFile, err := os.Open(txtPath)
	ErrCheck(err, "Cannot open input file: "+txtPath)
	defer txtFile.Close()

	// Scan rows of input txt
	var lines []string
	txtScanner := bufio.NewScanner(txtFile)
	lineNum := 0
	for txtScanner.Scan() {
		if lineNum == 0 {
			// Get all column names and deliminate by tab
			colNames := strings.Split(txtScanner.Text(), "\t")
			// Change column names in "bed" and "bed2"
			colNames = ChangeColName(colNames, renamemap)

			// Get column names until "Otherinfo1"
			var trimedCols []string
			for _, colName := range colNames {
				trimedCols = append(trimedCols, colName)
				if colName == "Otherinfo1" {
					break
				}
			}

			// Create a silice including VCF like columns
			vcflikeCols := []string{
				"Otherinfo2", "Otherinfo3", "CHROM", "POS", "ID",
				"REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
			}
			// Create a slice for header including all column names
			newCols := append(append(trimedCols, vcflikeCols...), samples...)
			// Delimitate by tab and add to lines
			lines = append(lines, strings.Join(newCols, "\t"))

		} else {
			// Non header lines are added to lines without any change
			lines = append(lines, txtScanner.Text())
		}
		lineNum++
	}

	if err := txtScanner.Err(); err != nil {
		panic(err)
	}
	return lines
}

// Get sample names from VCF
func GetSamples(vcfPath string) []string {
	vcfFile, err := os.Open(vcfPath)
	ErrCheck(err, "Cannot open input file: "+vcfPath)
	defer vcfFile.Close()

	// Scan rows of input VCF
	var samples []string
	vcfScanner := bufio.NewScanner(vcfFile)
	for vcfScanner.Scan() {
		line := vcfScanner.Text()
		// Search the row starting with '#CHROM'
		if strings.HasPrefix(line, "#CHROM") {
			columns := strings.Split(line, "\t")
			if len(columns) >= 10 {
				// Get sample names form samples info columns
				samples = columns[9:]
			} else {
				fmt.Println("This VCF does not have Sample information")
				os.Exit(1)
			}
			break
		}
	}
	if err := vcfScanner.Err(); err != nil {
		panic(err)
	}
	return samples
}

func SaveRenamedFile(lines []string, outPath string) {
	outputFile, err := os.Create(outPath)
	ErrCheck(err, "Cannot create output file: "+outPath)
	defer outputFile.Close()

	outputWriter := bufio.NewWriter(outputFile)
	for _, line := range lines {
		_, err := outputWriter.WriteString(line + "\n")
		if err != nil {
			panic(err)
		}
	}
	outputWriter.Flush()
}

func ErrCheck(err error, msg string) {
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s: %v\n", msg, err)
		os.Exit(1)
	}
}
