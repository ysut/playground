package main

import (
	// "bufio"
	"database/sql"
	// "fmt"
	"log"
	"os"
	// "strings"

	_ "github.com/mattn/go-sqlite3"
)

func main() {
	if len(os.Args) < 2 {
		log.Fatal("Usage: importtxt <dbfile> <txtfile>")
	}

	dbPath := os.Args[1]
	// txtFilePath := os.Args[2]

	db, err := sql.Open("sqlite3", dbPath)
	if err != nil {
		log.Fatalf("Error opening database: %v", err)
	}

	defer db.Close()

	// Create table
	createTableSQL := `CREATE TABLE IF NOT EXISTS my_table (
		id INTEGER PRIMARY KEY AUTOINCREMENT,
		column1 TEXT,
		column2 TEXT
	);`
	_, err = db.Exec(createTableSQL)
	if err != nil {
		log.Fatalf("Error creating table: %v", err)
	}


	
}