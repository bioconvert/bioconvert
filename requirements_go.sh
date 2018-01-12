go get -u github.com/golang/dep/cmd/dep

go get github.com/fredericlemoine/goalign/
cd $GOPATH/src/github.com/fredericlemoine/goalign/
dep ensure
make
cd -

go get github.com/fredericlemoine/gotree/
cd $GOPATH/src/github.com/fredericlemoine/gotree/
dep ensure
make
cd -

