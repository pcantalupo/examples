Reading large tables into R
===========================

Reading large tables from text files into R is possible but knowing a few tricks will make your life a lot easier and make R run a lot faster.

First, type `?read.table` in R and hit enter to read the help page for `read.table`. It contains many hints for how to read in large tables. Of course, help pages tend to be a little confusing so here is a distillation of a few of the golden rules.

The following options to `read.table()` can affect R's ability to read large tables:

- nrows
    - R doesn't know how many rows it's going to read in so it first makes a guess, and then when it runs out of room it allocates more memory. The constant allocations can take a lot of time, and if R overestimates the amount of memory it needs, your computer might run out of memory.
    - Of course, you may not know how many rows your table has. The easiest way to find this out is to use the `wc` command. So if you run `wc -l data.txt` on the command line, then it will report to you the number of lines in the file. You can then pass this number to the `nrows` argument of `read.table()` like this `nrows=231,238,977`.
    - If you can't use `wc` for some reason, but you know that there are definitely less than, say, `N` rows, then you can specify `nrows = N` and things will still be okay. A mild overestimate for `nrows` is better than none at all.
- colClasses
    - This option takes a vector whose length is equal to the number of columns in year table. Specifying this option instead of using the default can make `read.table` run MUCH faster, often twice as fast. In order to use this option, you have to know the class of each column in your data frame.
        - If all of the columns are "numeric", for example, then you can just set `colClasses = "numeric"`.
        - If you know that the first column is an integer, 2nd column is to be discarded and third is a factor then use `colClasses=c("integer","NULL","factor")`. If you don't need a column, specifying `NULL` will save lots of memory and time.
        - If you don't know classes of the columns, then you can have R do some of the work for you! Using the code below, I read in 10 rows of the table and determine the column class with `class` and then use the `classes` vector in `read.table`.
        ```
        tab5rows <- read.table("data.txt", header = TRUE, nrows = 10)
        classes <- sapply(tab5rows, class)
        tabAll <- read.table("data.txt", header = TRUE, colClasses = classes)
        ```
- commentChar
    - If your file has no comments in it (e.g. lines starting with `#`), then setting `comment.char = ""` will sometimes make `read.table()` run faster.

Other things to keep in mind:

- Use the save function to save intermediate results in `.RData` files. Compared to reading text files using `read.table`, `.RData` files are orders of magnitude faster to load back into R using the `load` function. Only save the things you need for computations later.
- Finally, avoid doing large vector operations when possible. For example, doing a bunch of things (calculate annotation error, segmentation models, etc) for each of 575 copy number profiles independently. If you store all the profiles in 1 data.frame with 4,616,846 rows and a column for `profile.id`, then calculations are much slower than if you split it into a list of 575 data.frames, each with 4000 rows on average. However the speedups using this trick are very much dependent on the number of rows in the big table and the number of groups you are splitting into, so you kind of have to try a few approaches and see what works fastest.

---
References:
* http://cbio.ensmp.fr/~thocking/reading-large-text-files-into-R.html
* http://www.biostat.jhsph.edu/~rpeng/docs/R-large-tables.html
