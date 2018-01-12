go get -u github.com/golang/dep/cmd/dep
go get github.com/fredericlemoine/goalign/
cd $GOPATH/src/github.com/fredericlemoine/goalign/
$GOPATH/bin/dep ensure
make
