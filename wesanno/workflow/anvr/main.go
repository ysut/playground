package main

import (
	"bufio"
	"fmt"
	"log"
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
		fmt.Println("Usage: anvrenamefilter <input.txt> <config.toml>")
		os.Exit(1)
	}

	// Read a config file
	configPath := os.Args[2]
	var config Config
	if _, err := toml.DecodeFile(configPath, &config); err != nil {
		fmt.Fprintf(os.Stderr, "Error reading config file: %v\n", err)
		os.Exit(1)
	}

	// #1. Get sample names from VCF
	txtPath := os.Args[1]
	vcfPath := strings.TrimSuffix(txtPath, ".txt") + ".vcf"
	var samples []string = GetSamples(vcfPath)

	// #2. Rename columns and extract exome variants
	var lines = RenameFilterProcess(txtPath, samples, config.RenameMap)

	// #3. Write a new txt
	outPath := strings.TrimSuffix(txtPath, ".txt") + ".renamed.txt"
	SaveRenamedFile(lines, outPath)
}

// Functions
// Have not been refactored yet -> ToDo: Refactoring these functions
func ChangeColName(colNames []string, renameMap map[string]string) []string {
	for i, colName := range colNames {
		for key, value := range renameMap {
			if colName == key {
				colNames[i] = value
			}
		}
	}
	return colNames
}

func GetColumnIndex(colNames []string, target string) int {
	for i, colName := range colNames {
		if colName == target {
			return i
		}
	}
	return -1
}

func GenerateNewColNames(colNames []string, samples []string) []string {
	var trimedCols []string
	for _, colName := range colNames {
		trimedCols = append(trimedCols, colName)
		if colName == "Otherinfo1" {
			break
		}
	}
	vcflikeCols := []string{
		"Otherinfo2", "Otherinfo3", "CHROM", "POS", "ID",
		"REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
	}
	newCols := append(append(trimedCols, vcflikeCols...), samples...)

	return newCols
}

func RenameFilterProcess(txtPath string, samples []string, renamemap map[string]string) []string {
	// Open a file annotated by ANNOVAR
	txtFile, err := os.Open(txtPath)
	ErrCheck(err, "Cannot open input file: "+txtPath)
	defer txtFile.Close()

	// Define a variable to store lines of the new txt
	var lines []string

	// Define variables to store the index of "CHROM" and "Func.refGene"
	// Intialize with -1
	var chromIndex, funcIndex int = -1, -1

	// Define slice to store the chromosomal values and the functional annotations
	chromValues := []string{
		"1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
		"11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
		"21", "22", "X", "Y",
	}

	functionalAnnotations := []string{
		"exonic", "splicing", "exonic;splicing",
	}

	// Scan rows of input txt
	txtScanner := bufio.NewScanner(txtFile)
	lineNum := 0

	for txtScanner.Scan() {
		var line string = txtScanner.Text()
		var rows []string = strings.Split(line, "\t")

		if lineNum == 0 {
			// Column names processing
			// 1st process: Change column names
			// Change column names as defined in the config file
			var colNames []string = ChangeColName(rows, renamemap)

			// 2nd process: Get the index of "CHROM" and "Func.refGene"
			chromIndex = GetColumnIndex(colNames, "Chr")
			funcIndex = GetColumnIndex(colNames, "Func.refGene")
			if chromIndex == -1 || funcIndex == -1 {
				log.Fatal("Cannot find 'CHROM' or 'Func.refGene' in the input file")
			}

			// 3rd process: Generate new column names
			var newCols []string = GenerateNewColNames(colNames, samples)

			// 4th process: Add new column names to lines
			lines = append(lines, strings.Join(newCols, "\t"))

		} else {
			// Variant rows processing
			chromValue := rows[chromIndex]
			funcValue := rows[funcIndex]

			// Pickup mitochondrial variants
			if chromValue == "MT" {
				lines = append(lines, line)
				continue
			}

			// Liner search for the chromValue and funcValue
			for _, funnanno := range functionalAnnotations {
				if funcValue == funnanno {
					for _, chrom := range chromValues {
						if chromValue == chrom {
							lines = append(lines, line)
							break
						}
					}
				}
			}
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
