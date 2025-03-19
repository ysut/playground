package main

import "fmt"

func Hello(name string) string {
	message := fmt.Sprintf("Hello, %s!", name)
	return message
}

func main() {
	message := Hello("World")
	fmt.Println(message)
}