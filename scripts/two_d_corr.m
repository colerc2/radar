function [ least_squares ] = two_d_corr( template, database )

number_of_comparisons = size(database,1) - size(template,1) + 1;
for ii = 1:number_of_comparisons
   diff = template - database(ii:ii+size(template,1)-1,:);
   least_squares(ii) = sum(sum(diff.^2));
end


end

