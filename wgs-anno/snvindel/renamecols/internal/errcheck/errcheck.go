package errcheck

import (
	"fmt"
	"os"
)

func ErrCheck(err error, msg string) {
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s: %v\n", msg, err)
		os.Exit(1)
	}
}
