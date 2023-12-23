function output = plotting_RDM(rdm, my_title)
output = 0;
figure
colormap(jet); % Or any other colormap of your choice
imagesc(rdm);
title(my_title)
xlabel('doppler bins');
ylabel('range bins');
colorbar; % Add a colorbar to the plot
end