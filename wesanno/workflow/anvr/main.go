package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"
	"github.com/BurntSushi/toml"
)

type Config struct {
    Columns map[string]string `toml:"columns"`
}


func main() {
	// Check number of args
	if len(os.Args) < 2 {
		fmt.Println("Please enter the input file path.")
		os.Exit(1)
	}
	configPath := "rename.toml"

	var config Config
	if _, err := toml.DecodeFile(configPath, &config); err != nil {
		fmt.Fprintf(os.Stderr, "Error reading config file: %v\n", err)
		os.Exit(1)
	}

	// #1. Get sample names from VCF
	vcfPath := strings.TrimSuffix(os.Args[1], ".txt") + ".vcf"
	var samples []string = GetSamples(vcfPath)

	fmt.Println(samples)

	// for key, value := range config.Columns {
	// 	fmt.Printf("Column key: %s, Column name: %s\n", key, value)



	// }
	fmt.Println(config.Columns["bed"])



	// // #2. Rename columns in txt
	// // Rename columns and return lines to write a new txt
	// txtPath := os.Args[1]
	// var lines []string = ProcessTxtFile(txtPath, samples, bedColName, bed2ColName)

	// // #3. Write a new txt
	// outPath := strings.TrimSuffix(txtPath, ".txt") + ".renamed.txt"
	// SaveRenamedFile(lines, outPath)
}		


func ChangeBedColName(colNames []string, bed string, bed2 string) []string {
	for i, colName := range colNames {
		if colName == "bed" {
			colNames[i] = bed
		} else if colName == "bed2" {
			colNames[i] = bed2
		}
	}
	return colNames
}

func ProcessTxtFile(txtPath string, samples []string, bed string, bed2 string) []string {
	// Open genome_summary.txt
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
			colNames = ChangeBedColName(colNames, bed, bed2)

			// Get column names until "Otherinfo1"
			var trimedCols []string
			for _, colName := range colNames {
				trimedCols = append(trimedCols, colName)
				if colName == "Otherinfo1" {
					break
				}
			}

			// Create a silice including VCF like columns
			vcfSampleCols := []string{
				"Otherinfo2", "Otherinfo3", "CHROM", "POS", "ID",
				"REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
			}

			// Create a slice for header including all column names
			newCols := append(append(trimedCols, vcfSampleCols...), samples...)
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







func ErrCheck(err error, msg string) {
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s: %v\n", msg, err)
		os.Exit(1)
	}
}