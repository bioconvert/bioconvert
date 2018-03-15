[ -z "$TRAVIS_PYTHON_VERSION" ] || GOROOT=$(dirname $(which conda))/../go && unset GOPATH
[ -z "$GOPATH" ] && GOPATH="$HOME/go/"
PATH=$GOPATH/bin:$GOROOT:$PATH
go get -u github.com/golang/dep/cmd/dep
go get -u github.com/fredericlemoine/gotree/
cd $GOPATH/src/github.com/fredericlemoine/gotree/
dep ensure
make
