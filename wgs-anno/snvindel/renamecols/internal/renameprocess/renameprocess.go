package renameprocess

import (
	"bufio"
	"fmt"
	"os"
	"renamecols/internal/errcheck"
	"strings"
)

// Change column names in "bed" and "bed2"
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

// Get sample names from VCF
func GetSamples(vcfPath string) []string {
	vcfFile, err := os.Open(vcfPath)
	errcheck.ErrCheck(err, "Cannot open input file: "+vcfPath)
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

func ProcessTxtFile(txtPath string, samples []string, bed string, bed2 string) []string {
	// Open genome_summary.txt
	txtFile, err := os.Open(txtPath)
	errcheck.ErrCheck(err, "Cannot open input file: "+txtPath)
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

func SaveRenamedFile(lines []string, outPath string) {
	outputFile, err := os.Create(outPath)
	errcheck.ErrCheck(err, "Cannot create output file: "+outPath)
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
