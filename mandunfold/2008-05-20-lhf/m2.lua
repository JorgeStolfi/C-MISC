local N=1024
print("P2\n",N,N,255)
for i=1,N do
	local x=-2+4*(i-1)/N
	for j=1,N do
		local y=-2+4*(j-1)/N
		--local xc,yc = x^2-y^2, 2*x*y
		--local xc,yc = x,y
		local xc,yc = 0.25-(x^2-y^2), -2*x*y
		local X,Y=0,0
		local p=0
		for n=1,255 do
			X,Y = X^2-Y^2+xc, 2*X*Y+yc
			if X^2+Y^2>4 then p=n break end
		end
		print(p)
	end
end
