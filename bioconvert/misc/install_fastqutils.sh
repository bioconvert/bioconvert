[ -z "$TRAVIS_PYTHON_VERSION" ] || GOROOT=$(dirname $(which conda))/../go && unset GOPATH
[ -z "$GOPATH" ] && GOPATH="$HOME/go/"
PATH=$GOPATH/bin:$GOROOT:$PATH
go get -u github.com/fredericlemoine/fastqutils
cd $GOPATH/src/github.com/fredericlemoine/fastqutils/
make
