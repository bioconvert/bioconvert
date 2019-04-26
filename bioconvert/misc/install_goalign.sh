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
[ -z "$TRAVIS_PYTHON_VERSION" ] || GOROOT="$HOME/miniconda3/envs/testenv/go" && unset GOPATH
[ -z "$GOPATH" ] && GOPATH="$HOME/go/"
PATH=$GOPATH/bin:$GOROOT:$PATH
go get -u github.com/golang/dep/cmd/dep
go get -u github.com/evolbioinfo/goalign/
cd $GOPATH/src/github.com/evolbioinfo/goalign/
$GOPATH/bin/dep ensure
# the goalign Mafile consider GOPATH to be empty on some platform
# so let us provide it when calling make

make GOPATH="$GOPATH" && make install
