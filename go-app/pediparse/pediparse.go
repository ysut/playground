
package main

import (
    "bufio"
    "fmt"
    "os"
    "strings"
    "strconv"
)

func main() {
    scanner := bufio.NewScanner(os.Stdin)

    for scanner.Scan() {
        line := scanner.Text()
        if strings.TrimSpace(line) == "" {
            continue
        }
        cols := strings.Split(line, "\t")
        if len(cols) < 6 {
            // PED format 6列以上必須
            continue
        }
        family := cols[0]
        indiv  := cols[1]
        father := cols[2]
        mother := cols[3]
        sex    := cols[4]
        pheno  := cols[5]

        // ここで続柄を判定するロジックを入れる
        // 例: father != "0" or mother != "0" → child?
        //     sex = 1 & 誰かのfather → father?
        //     独自ルール: IDが若い方をproband etc.
        // 以下はダミーで一律 "role=?" を付ける
        role := "?"

        // 例: fatherId=0, motherId=0なら founder
        if father == "0" && mother == "0" {
            role = "founder"
        }
        // 例: phenotype=2 の人は proband とする（ダミーの例）
        phenoInt, _ := strconv.Atoi(pheno)
        if phenoInt == 2 {
            role = "proband"
        }

        // タブ区切りで出力
        // 例: 既存列 + role列
        fmt.Printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
            family, indiv, father, mother, sex, pheno, role)
    }

    if err := scanner.Err(); err != nil {
        fmt.Fprintf(os.Stderr, "Error: %v\n", err)
        os.Exit(1)
    }
}
