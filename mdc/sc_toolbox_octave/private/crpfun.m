function f = crpfun(x,fdat)
%CRPFUN (not intended for calling directly by the user)
%   Nonlinear function for CRPARAM.
%   Copyright 1999 by Toby Driscoll.
%   $Id: crpfun.m 298 2009-09-15 14:36:37Z driscoll $

[n,beta,crtarget,Q,qdat,nproc] = deal(fdat{:});

crprever = exp(x);			% prevertex crossratios

use_parallel = nproc > 0 && exist('parcellfun') == 2;

if use_parallel
	crimage = cell2mat(parcellfun(nproc,@(k) crpfun_one(k,crprever,Q,beta,qdat), ...
		num2cell(1:n-3),'UniformOutput',false,'VerboseLevel',1, ...
		'ErrorHandler',@dbg_err_handler));
	crimage = crimage(:);
else
	crimage = zeros(n-3,1);			% image vertex crossratios
	% Compute crossratio for each image quadrilateral
	for k = 1:n-3
		crimage(k) = crpfun_one(k,crprever,Q,beta,qdat);
	end
end

% Logarithmic scaling for residual
f = log(abs(crimage./crtarget));
end

function val = dbg_err_handler(err,varargin)
	fprintf(2,'WORKER ERROR: %s\n',err.message);
	val = NaN;
end

function val = crpfun_one(k,crprever,Q,beta,qdat)
	prever = crembed(crprever,Q,k);
	w = -crquad(prever(Q.qlvert(:,k)),Q.qlvert(:,k),prever,beta,qdat);
	val = (w(2)-w(1))*(w(4)-w(3))/((w(3)-w(2))*(w(1)-w(4)));
end
