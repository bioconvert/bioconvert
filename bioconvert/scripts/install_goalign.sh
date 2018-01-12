[ -z "$GOPATH" ] && GOPATH="$HOME/go/"
PATH=$GOPATH/bin:$PATH
go get -u github.com/golang/dep/cmd/dep
go get github.com/fredericlemoine/goalign/
cd $GOPATH/src/github.com/fredericlemoine/goalign/
dep ensure
make
