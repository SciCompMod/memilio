import cv2
import numpy as np
import pandas as pd


def impute_comix(matrix, env, path):
    polymod = pd.read_csv(path + 'Polymod' + env + '.csv', index_col=0).values

    eig_polymod, _ = np.linalg.eig(polymod[2:, 2:])
    eig_comix, _ = np.linalg.eig(matrix[2:, 2:])

    max_eig_polymod = np.max(eig_polymod)
    max_eig_comix = np.max(eig_comix)

    ratio = max_eig_comix / max_eig_polymod

    new_comix = np.zeros((8, 8))

    for i in range(8):
        for j in range(2):
            new_comix[i, j] = polymod[i, j] * ratio

    new_comix[:,2:] = matrix[:,2:]

    return new_comix

def read_from_screenshot(image, env, path, impute = False, show_result = False):
    print(path + image + env + '.png')
    image_values = np.flip(cv2.imread(path + image + env + '.png'), axis=0)
    colorbar = np.flip(cv2.imread(path + 'Colorbar.png'), axis=0)




    colors = np.zeros((8,8,3))

    for i in range(8):
        for j in range(8):
            colors[i,j,:] = image_values[int(image_values.shape[0]*i/8 + image_values.shape[0]/16),
                            int(image_values.shape[1]*j/8 + image_values.shape[1]/16), :]

    print(colorbar.shape)


    contact_matrix = np.zeros((8,8))
    for l in range(8):
        for k in range(8):
            temp = np.abs(colorbar[:,int(colorbar.shape[1]/2),:] - colors[l,k,:])

            color_norm = np.zeros(colorbar.shape[0])
            for i in range(len(colorbar[:,0,0])):
                color_norm[i] = np.linalg.norm(temp[i,:])
            contact_matrix[l,k] = np.argmin(color_norm)*9/colorbar.shape[0]

    if impute:
        contact_matrix = impute_comix(contact_matrix, env, path)

    print(contact_matrix)
    columns = ['0-4', '5-17', '18-29', '30-39', '40-49', '50-59', '60-69', '70+']
    df = pd.DataFrame(data = contact_matrix, index=columns, columns=columns)
    df.to_csv(path + image + env + '.csv')

    if show_result:
        cv2.imshow('image', image_values)
        scale_percent = 10000  # percent of original size
        width = int(contact_matrix.shape[1] * scale_percent / 100)
        height = int(contact_matrix.shape[0] * scale_percent / 100)
        dim = (width, height)
        # resize image
        resized = cv2.resize(contact_matrix, dim, interpolation=cv2.INTER_AREA)

        print('Resized Dimensions : ', resized.shape)

        cv2.imshow("Resized image", resized/10)
        cv2.waitKey(0)
        cv2.destroyAllWindows()


if __name__ == "__main__":
    path = '../../../data/contact_data/images/'
    read_from_screenshot('Polymod', '', path, show_result=False, impute=False)
    read_from_screenshot('Comix', '', path, show_result=False, impute=False)
    names = ['_Physical', '_Home', '_Work', '_Other', '_School']
    for env in names:
        read_from_screenshot('Polymod', env, path, show_result=False, impute=False)
        read_from_screenshot('Comix', env, path, show_result=False, impute=True)
