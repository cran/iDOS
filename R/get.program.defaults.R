get.program.defaults <- function() {

	# make a list of potential locations for the datasets file
	program.data.dirs <- paste(system.file("programdata/", package = "iDOS"), "/", sep = "");
	
	return (
		list(
			"program.data.dir" = program.data.dirs[1],
			"test.data.dir" = paste(program.data.dirs[1], "testdata/", sep = "")
			)
		);
	}
