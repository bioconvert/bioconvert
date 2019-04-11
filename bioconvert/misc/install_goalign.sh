###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Authors: see CONTRIBUTORS.rst                                           #
# Copyright © 2018  Institut Pasteur, Paris and CNRS.                     #
# See the COPYRIGHT file for details                                      #
#                                                                         #
# bioconvert is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# bioconvert is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License       #
# along with this program (COPYING file).                                 #
# If not, see <http://www.gnu.org/licenses/>.                             #
###########################################################################
# GOPATH workspace directory
# A priori no need to set GOROOT but need to be set if binary distribution
# is not in the default location.
# GOPATH can contains several entries like PATH
echo $HOME

# $(which conda) may not give the correct path (april 2019) replaced by harcoded
# value miniconda3/testenv/go since a deidcated environement is set
#[ -z "$TRAVIS_PYTHON_VERSION" ] || GOROOT="$(dirname $(which conda))"/../go && unset GOPATH
[ -z "$TRAVIS_PYTHON_VERSION" ] || GOROOT="$HOME/miniconda3/testenv/go" && unset GOPATH
[ -z "$GOPATH" ] && GOPATH="$HOME/go/"

echo $GOPATH
echo $GOROOT

PATH=$GOPATH/bin:$GOROOT:$PATH

echo $PATH

go get -u github.com/golang/dep/cmd/dep
go get -u github.com/evolbioinfo/goalign/
cd $GOPATH/src/github.com/evolbioinfo/goalign/
$GOPATH/bin/dep ensure
# the goalign Mafile consider GOPATH to be empty on some platform
# so let us provide it when calling make
make GOPATH="$GOPATH"
