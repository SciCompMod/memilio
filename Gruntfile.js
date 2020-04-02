
module.exports = function(grunt) {
    // load grunt tasks based on dependencies in package.json
    require('load-grunt-tasks')(grunt);

    grunt.config.init({
      useminPrepare: {
          html: 'index.html',
          options: {
            dest: 'dist'
          }
      },
      copy:{
        html: {
            src: './index.html', dest: 'dist/index.html',
        },
      main: {
        files: [
          {expand: true, src: ['./*.csv'], dest: './dist/', filter: 'isFile'},
          {expand: true, src: ['./fonts/*'], dest: './dist/', filter: 'isFile'},
          {expand: true, src: ['./images/*'], dest: './dist/', filter: 'isFile'},
          {expand: true, src: ['./data/*'], dest: './dist/', filter: 'isFile'},
          {expand: true, src: ['./js/*.js'], dest: './dist/'},
          {expand: true, src: ['./*.js'], dest: './dist/'},
          {expand: true, src: ['./*.css'], dest: './dist/'},
        ]
      }
      }
    });

    grunt.registerTask('default',[
        'copy:html',
        'copy:main',
    ]);
}
